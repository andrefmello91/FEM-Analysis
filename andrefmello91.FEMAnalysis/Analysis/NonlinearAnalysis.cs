using System;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Nonlinear analysis class
	/// </summary>
	public class NonlinearAnalysis : Analysis<IFiniteElement>
	{

		#region Fields

		/// <summary>
		///     The list of iteration results.
		/// </summary>
		/// <remarks>
		///     This is cleaned at the beginning of a load step.
		/// </remarks>
		private readonly List<IterationResult> _iterations = new();

		/// <summary>
		///     The list of load step results.
		/// </summary>
		private readonly List<LoadStepResult> _loadSteps = new();

		/// <summary>
		///     Field to store the DoF index for monitored displacements.
		/// </summary>
		private int? _monitoredIndex;

		#endregion

		#region Properties

		/// <inheritdoc />
		/// <remarks>
		///     The displacements of current iteration.
		/// </remarks>
		public override Vector<double>? DisplacementVector
		{
			get => OngoingIteration.Displacements;
			protected set
			{
				if (value is null)
					return;

				OngoingIteration.Displacements = value;
			}
		}

		/// <inheritdoc />
		/// <remarks>
		///     The stiffness of current iteration.
		/// </remarks>
		public override Matrix<double>? GlobalStiffness
		{
			get => OngoingIteration.Stiffness;
			protected set
			{
				if (value is null)
					return;

				OngoingIteration.Stiffness = value;
			}
		}

		/// <summary>
		///     Get/set the maximum number of iterations.
		/// </summary>
		public int MaxIterations { get; set; }

		/// <summary>
		///     Get/set the minimum number of iterations.
		/// </summary>
		public int MinIterations { get; set; }

		/// <summary>
		///     Get/set the number of load steps to execute.
		/// </summary>
		public int NumLoadSteps { get; set; }

		/// <summary>
		///     The nonlinear equation solver.
		/// </summary>
		public NonLinearSolver Solver { get; set; }

		/// <summary>
		///     Get/set when to stop analysis.
		/// </summary>
		public bool Stop { get; private set; }

		/// <summary>
		///     Get/set the stop message.
		/// </summary>
		public string StopMessage { get; private set; } = string.Empty;

		/// <summary>
		///     Get/set the convergence tolerance.
		/// </summary>
		public double Tolerance { get; set; }

		/// <summary>
		///     The current load step result.
		/// </summary>
		private LoadStepResult CurrentLoadStep => _loadSteps[^1];

		/// <summary>
		///     The results of the current solution (last solved iteration [i - 1]).
		/// </summary>
		private IterationResult CurrentSolution => _iterations[^2];

		/// <summary>
		///     The last load step result.
		/// </summary>
		private LoadStepResult LastLoadStep => _loadSteps.Count > 1
			? _loadSteps[^2]
			: CurrentLoadStep;

		/// <summary>
		///     The results of the last solution (penultimate solved iteration [i - 2]).
		/// </summary>
		private IterationResult LastSolution => _iterations[^3];

		/// <summary>
		///     Get the load factor of the current iteration.
		/// </summary>
		private double LoadFactor => (double) (int) CurrentLoadStep / NumLoadSteps;

		/// <summary>
		///     The results of the ongoing iteration.
		/// </summary>
		private IterationResult OngoingIteration => _iterations[^1];

		/// <summary>
		///     Get/set the internal force vector of current iteration.
		/// </summary>
		private Vector<double> InternalForces
		{
			get => OngoingIteration.InternalForces;
			set => OngoingIteration.UpdateForces(CurrentLoadStep.Forces, value);
		}

		/// <summary>
		///     Get the residual force vector of current iteration.
		/// </summary>
		private Vector<double> ResidualForces => OngoingIteration.ResidualForces;

		#endregion

		#region Constructors

		/// <summary>
		///     Nonlinear analysis constructor.
		/// </summary>
		/// <param name="nonlinearInput">The <see cref="IFEMInput{INonlinearElement}" />.</param>
		/// <param name="solver">The <see cref="NonLinearSolver" /> to use.</param>
		/// <param name="numLoadSteps">The number of load steps to perform (default: 50).</param>
		/// <param name="tolerance">The convergence tolerance (default: 1E-3).</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step (default: 10000).</param>
		/// <param name="minIterations">Minimum number of iterations for each load step (default: 2).</param>
		public NonlinearAnalysis(
			IFEMInput<IFiniteElement> nonlinearInput,
			NonLinearSolver solver = NonLinearSolver.NewtonRaphson,
			int numLoadSteps = 50,
			double tolerance = 1E-3,
			int maxIterations = 10000,
			int minIterations = 2)
			: base(nonlinearInput)
		{
			Solver        = solver;
			NumLoadSteps  = numLoadSteps;
			Tolerance     = tolerance;
			MaxIterations = maxIterations;
			MinIterations = minIterations;
		}

		#endregion

		#region Methods

		/// <summary>
		///     Calculate the secant stiffness increment.
		/// </summary>
		/// <param name="currentStiffness">The stiffness matrix from current iteration.</param>
		/// <param name="currentDisplacements">The displacement vector from current iteration.</param>
		/// <param name="lastDisplacements">The displacement vector from the last iteration.</param>
		/// <param name="currentResidual">The residual force vector from current iteration.</param>
		/// <param name="lastResidual">The residual force vector from last iteration.</param>
		/// <returns>
		///     <see cref="Matrix{T}" />
		/// </returns>
		public static Matrix<double> SecantIncrement(Matrix<double> currentStiffness, Vector<double> currentDisplacements, Vector<double> lastDisplacements, Vector<double> currentResidual, Vector<double> lastResidual)
		{
			// Calculate the variation of displacements and residual as vectors
			Vector<double>
				dDisp = currentDisplacements - lastDisplacements,
				dRes  = currentResidual - lastResidual;

			return
				((dRes - currentStiffness * dDisp) / dDisp.Norm(2)).ToColumnMatrix() * dDisp.ToRowMatrix();
		}

		/// <summary>
		///     Calculate the tangent stiffness increment.
		/// </summary>
		/// <param name="currentInternalForces">The internal force vector from current iteration.</param>
		/// <param name="lastInternalForces">The internal force vector from last iteration.</param>
		/// <param name="currentDisplacements">The displacement vector from current iteration.</param>
		/// <param name="lastDisplacements">The displacement vector from the last iteration.</param>
		/// <returns>
		///     <see cref="Matrix{T}" />
		/// </returns>
		public static Matrix<double> TangentIncrement(Vector<double> currentInternalForces, Vector<double> lastInternalForces, Vector<double> currentDisplacements, Vector<double> lastDisplacements)
		{
			// Get variations
			var dF = currentInternalForces - lastInternalForces;
			var dU = currentDisplacements  - lastDisplacements;
			
			return
				dF.ToColumnMatrix() * dU.ToRowMatrix();
		}
		
		/// <summary>
		///     Calculate the force based convergence.
		/// </summary>
		/// <param name="residualForces">The residual forces of the current iteration.</param>
		/// <param name="appliedForces">The applied forces of the current load step.</param>
		internal static double ForceConvergence(IEnumerable<double> residualForces, IEnumerable<double> appliedForces)
		{
			double
				num = residualForces.Sum(n => n * n),
				den = 1 + appliedForces.Sum(n => n * n);

			return
				num / den;
		}

		/// <summary>
		///     Execute the analysis.
		/// </summary>
		/// <inheritdoc cref="Initiate" />
		/// <inheritdoc cref="StepAnalysis" />
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis{TFiniteElement}.ForceVector" /> (default: 1).</param>
		public void Execute(int? monitoredIndex = null, bool simulate = false, double loadFactor = 1)
		{
			// Get force vector
			ForceVector = FemInput.ForceVector * loadFactor;

			// Initiate lists
			Initiate(monitoredIndex);

			// Analysis by load steps
			StepAnalysis(simulate);

			// Set Reactions
			FemInput.Grips.SetReactions(GetReactions());
		}

		/// <summary>
		///     Generate an <see cref="FEMOutput" /> from analysis results.
		/// </summary>
		/// <returns>
		///     null if no monitored index was provided.
		/// </returns>
		public FEMOutput GenerateOutput() => new(_loadSteps);

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix{T}" /> of current iteration.
		/// </summary>
		/// <param name="simplify">Simplify stiffness?</param>
		protected override void UpdateStiffness(bool simplify = true)
		{
			switch (Solver)
			{
				case NonLinearSolver.Secant:
					// Increment current stiffness
					GlobalStiffness += SecantIncrement(CurrentSolution.Stiffness, CurrentSolution.Displacements, LastSolution.Displacements, CurrentSolution.ResidualForces, LastSolution.ResidualForces);
					break;

				// For Newton-Raphson
				case NonLinearSolver.NewtonRaphson:
				case NonLinearSolver.ModifiedNewtonRaphson when (int) OngoingIteration == 1:
					// // Update stiffness in elements
					// FemInput.Elements.UpdateStiffness();
					//
					// // Set new values
					// GlobalStiffness = FemInput.AssembleStiffness();
					
					GlobalStiffness += TangentIncrement(CurrentSolution.InternalForces, LastSolution.InternalForces, CurrentSolution.Displacements, LastSolution.Displacements);

					break;
				
				default:
					return;
			}

			// Simplify
			if (simplify)
				Simplify(GlobalStiffness, null, FemInput.ConstraintIndex);
		}

		/// <summary>
		///     Correct results from last load step after not achieving convergence.
		/// </summary>
		private void CorrectResults()
		{
			// Set displacements from last (current now) load step
			DisplacementVector = LastLoadStep.Displacements;
			GlobalStiffness    = LastLoadStep.Stiffness;
			FemInput.Grips.SetDisplacements(DisplacementVector);
			FemInput.Elements.UpdateDisplacements();

			// Calculate element forces
			FemInput.Elements.CalculateForces();
		}

		/// <summary>
		///     Initiate fields.
		/// </summary>
		/// <param name="monitoredIndex">The index of a degree of freedom to monitor, if wanted.</param>
		private void Initiate(int? monitoredIndex)
		{
			_monitoredIndex = monitoredIndex;

			// Initiate lists solution values
			_loadSteps.Clear();
			_loadSteps.Add(new LoadStepResult(ForceVector / NumLoadSteps, 1));

			_iterations.Clear();
			for (var i = 0; i < 3; i++)
				_iterations.Add(new IterationResult(FemInput.NumberOfDoFs));

			// Get the initial stiffness and force vector simplified
			GlobalStiffness = FemInput.AssembleStiffness();
			Simplify(GlobalStiffness, ForceVector, FemInput.ConstraintIndex);

			// Calculate initial displacements
			var fi = CurrentLoadStep.Forces;
			DisplacementVector = GlobalStiffness!.Solve(fi);

			// Update displacements in grips and elements
			FemInput.Grips.SetDisplacements(DisplacementVector);
			FemInput.Elements.UpdateDisplacements();

			// Calculate element forces
			FemInput.Elements.CalculateForces();

			// Update internal forces
			InternalForces = FemInput.AssembleInternalForces();
		}

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		private void Iterate()
		{
			ClearIterations();
			
			// Initiate first iteration
			OngoingIteration.Number = 0;

			// Iterate
			do
			{
				// Add iteration
				_iterations.Add(OngoingIteration.Clone());

				// Increase iteration count
				OngoingIteration.Number++;

				// Update stiffness and displacements
				UpdateDisplacements();
				UpdateStiffness();

				if ((int) OngoingIteration == 50)
				{
					
				}
					// Calculate element forces
				FemInput.Elements.CalculateForces();

				// Update internal forces
				InternalForces = FemInput.AssembleInternalForces();

				// Calculate convergence
				OngoingIteration.CalculateConvergence(CurrentLoadStep.Forces);
				
			} while (!IterativeStop());
		}

		/// <summary>
		///		Clear the iterations lists.
		/// </summary>
		protected virtual void ClearIterations()
		{
			// Clear iteration list
			if ((int) CurrentLoadStep <= 1 || (int) OngoingIteration < 4)
				return;
			
			_iterations.RemoveRange(..^2);
		}
		
		/// <summary>
		///     Check if iterative procedure must stop by achieving convergence or achieving the maximum number of iterations.
		/// </summary>
		/// <remarks>
		///     If the maximum number of iterations is reached, <see cref="Stop" /> is set to true.
		/// </remarks>
		/// <returns>
		///     True if convergence is reached or the maximum number of iterations is reached.
		/// </returns>
		private bool IterativeStop()
		{
			// Check if one stop condition is reached
			Stop = (int) OngoingIteration >= MaxIterations || ResidualForces.ContainsNaNOrInfinity() ||
			       DisplacementVector!.ContainsNaNOrInfinity() || GlobalStiffness!.ContainsNaN();

			switch (Stop)
			{
				// Check if maximum number of iterations is reached
				case true:
					StopMessage = $"Convergence not reached at load step {(int) CurrentLoadStep}";
					return Stop;

				default:
					return
						VerifyConvergence(OngoingIteration.Convergence);
			}
		}

		/// <summary>
		///     Save load step results after achieving convergence.
		/// </summary>
		private void SaveLoadStepResults()
		{
			CurrentLoadStep.IsCalculated  = true;
			CurrentLoadStep.Convergence   = OngoingIteration.Convergence;
			CurrentLoadStep.Displacements = DisplacementVector!;
			CurrentLoadStep.Stiffness     = GlobalStiffness!;

			if (!_monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = Length.FromMillimeters(DisplacementVector![_monitoredIndex.Value]);

			// Set to load step
			CurrentLoadStep.MonitoredDisplacement = new MonitoredDisplacement(disp, LoadFactor);
		}

		/// <summary>
		///     Execute step by step analysis.
		/// </summary>
		/// <param name="simulate">Set true to execute analysis until convergence is not achieved (structural failure).</param>
		private void StepAnalysis(bool simulate)
		{
			// Initiate first load step
			CurrentLoadStep.Number = 1;

			// Step-by-step analysis
			do
			{
				// Get the force vector
				CurrentLoadStep.Forces = LoadFactor * ForceVector;

				// Iterate
				Iterate();

				// Verify if convergence was not reached
				if (Stop)
					goto CorrectResults;

				// Set load step results
				SaveLoadStepResults();

				// Create load step
				_loadSteps.Add(CurrentLoadStep.Clone());

				// Increment load step
				CurrentLoadStep.Number++;
			} while (simulate || (int) CurrentLoadStep <= NumLoadSteps);

			CorrectResults:
			CorrectResults();
		}

		/// <summary>
		///     Update displacements.
		/// </summary>
		private void UpdateDisplacements()
		{
			// Increment displacements
			DisplacementVector -= GlobalStiffness!.Solve(CurrentSolution.ResidualForces);

			// Update displacements in grips and elements
			FemInput.Grips.SetDisplacements(DisplacementVector);
			FemInput.Elements.UpdateDisplacements();
		}

		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="convergence">
		///     Calculated convergence.
		///     <para>See: <see cref="ForceConvergence" />.</para>
		/// </param>
		private bool VerifyConvergence(double convergence) => convergence <= Tolerance && (int) OngoingIteration >= MinIterations;

		#endregion

	}
}