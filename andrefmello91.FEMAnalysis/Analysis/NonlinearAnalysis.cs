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
		///     This is cleaned at the beginning of a step.
		/// </remarks>
		private readonly List<IterationResult> _iterations = new();

		/// <summary>
		///     The list of step results.
		/// </summary>
		private readonly List<StepResult> _steps = new();

		/// <summary>
		///     Field to store the DoF index for monitored displacements.
		/// </summary>
		private int? _monitoredIndex;

		/// <summary>
		///		Get the displacement increment for the first iteration of the current step.
		/// </summary>
		private Vector<double> FirstIncrement => _iterations.Find(i => (int) i == 1)!.DisplacementIncrement;

		/// <summary>
		///		The initial parameter for calculating the current stiffness parameter.
		/// </summary>
		private double _k0;
		
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
		/// <remarks>
		///		Default: 10000
		/// </remarks>
		public int MaxIterations { get; set; } = 10000;

		/// <summary>
		///     Get/set the minimum number of iterations.
		/// </summary>
		/// <remarks>
		///		Default: 2
		/// </remarks>
		public int MinIterations { get; set; } = 2;

		/// <summary>
		///     Get/set the number of steps to execute.
		/// </summary>
		/// <remarks>
		///		Default: 50
		/// </remarks>
		public int NumberOfSteps { get; set; } = 50;

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
		///     Get/set the convergence tolerance for residual forces.
		/// </summary>
		/// <remarks>
		///		Default: 1E-3
		/// </remarks>
		public double ForceTolerance { get; set; } = 1E-3;
		
		/// <summary>
		///     Get/set the convergence tolerance for displacement increments.
		/// </summary>
		/// <remarks>
		///		Default: 1E-8
		/// </remarks>
		public double DisplacementTolerance { get; set; } = 1E-8;

		/// <summary>
		///     The current step result.
		/// </summary>
		private StepResult CurrentStep => _steps[^1];

		/// <summary>
		///     The results of the current solution (last solved iteration [i - 1]).
		/// </summary>
		private IterationResult CurrentSolution => _iterations[^2];

		/// <summary>
		///     The last step result.
		/// </summary>
		private StepResult LastStep => _steps.Count > 1
			? _steps[^2]
			: CurrentStep;

		/// <summary>
		///     The results of the last solution (penultimate solved iteration [i - 2]).
		/// </summary>
		private IterationResult LastSolution => _iterations[^3];

		/// <summary>
		///     Get the load factor of the current step.
		/// </summary>
		private double LoadFactor { get; set; }

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
			set => OngoingIteration.UpdateForces(SimplifiedForces(CurrentStep.Forces, FemInput.ConstraintIndex), value);
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
		public NonlinearAnalysis(IFEMInput<IFiniteElement> nonlinearInput, NonLinearSolver solver = NonLinearSolver.NewtonRaphson)
			: base(nonlinearInput)
		{
			Solver = solver;
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
				dU = currentDisplacements - lastDisplacements,
				dR = currentResidual - lastResidual;

			return
				((dR - currentStiffness * dU) / dU.Norm(2)).ToColumnMatrix() * dU.ToRowMatrix();
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
		///     Calculate the convergence.
		/// </summary>
		/// <param name="numerator">The residual forces of the current iteration or the initial displacement increment of the current step.</param>
		/// <param name="denominator">The applied forces of the current step or the displacement increment of the current iteration.</param>
		internal static double CalculateConvergence(IEnumerable<double> numerator, IEnumerable<double> denominator)
		{
			double
				num = numerator.Sum(n => n * n),
				den = 1 + denominator.Sum(n => n * n);

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

			// Analysis by steps
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
		public FEMOutput GenerateOutput() => new(_steps);

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix{T}" /> of current iteration.
		/// </summary>
		protected override void UpdateStiffness()
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
					// Update stiffness in elements
					FemInput.UpdateStiffness();
					
					// Set new values
					GlobalStiffness = FemInput.AssembleStiffness();
					
					// GlobalStiffness += TangentIncrement(CurrentSolution.InternalForces, LastSolution.InternalForces, CurrentSolution.Displacements, LastSolution.Displacements);

					break;
				
				default:
					return;
			}
		}

		/// <summary>
		///     Correct results from last step after not achieving convergence.
		/// </summary>
		private void CorrectResults()
		{
			// Set displacements from last (current now) step
			DisplacementVector = LastStep.Displacements;
			GlobalStiffness    = LastStep.Stiffness;
			FemInput.Grips.SetDisplacements(DisplacementVector);
			FemInput.UpdateDisplacements();

			// Calculate element forces
			FemInput.CalculateForces();
		}

		/// <summary>
		///     Initiate fields.
		/// </summary>
		/// <param name="monitoredIndex">The index of a degree of freedom to monitor, if wanted.</param>
		private void Initiate(int? monitoredIndex)
		{
			_monitoredIndex = monitoredIndex;
			LoadFactor      = StepIncrement(false);
			
			// Initiate lists solution values
			_steps.Clear();
			_steps.Add(new StepResult(LoadFactor * ForceVector, 1));

			_iterations.Clear();
			for (var i = 0; i < 3; i++)
				_iterations.Add(new IterationResult(FemInput.NumberOfDoFs));

			// Get the initial stiffness and force vector simplified
			GlobalStiffness = FemInput.AssembleStiffness();
			var stiffness = SimplifiedStiffness(GlobalStiffness, FemInput.ConstraintIndex);

			// Calculate initial displacements
			var fi = SimplifiedForces(CurrentStep.Forces, FemInput.ConstraintIndex);
			DisplacementVector = stiffness.Solve(fi);

			// Update displacements in grips and elements
			FemInput.Grips.SetDisplacements(DisplacementVector);
			FemInput.UpdateDisplacements();

			// Calculate element forces
			FemInput.CalculateForces();

			// Update internal forces
			InternalForces = FemInput.AssembleInternalForces();
		}

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		/// <inheritdoc cref="StepAnalysis"/>
		private void Iterate(bool simulate)
		{
			ClearIterations();

			// Add iteration
			_iterations.Add(OngoingIteration.Clone());

			// Initiate first iteration
			OngoingIteration.Number = 0;

			// Do initial steps
			PredictorStep(simulate);
			
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
				FemInput.CalculateForces();

				// Update internal forces
				InternalForces = FemInput.AssembleInternalForces();

				// Calculate convergence
				OngoingIteration.CalculateForceConvergence(CurrentStep.Forces);
				OngoingIteration.CalculateDisplacementConvergence(FirstIncrement);

			} while (!IterativeStop());
		}

		/// <summary>
		///		Steps to perform at the beginning of a load step.
		/// </summary>
		/// <inheritdoc cref="StepAnalysis"/>
		private void PredictorStep(bool simulate)
		{
			// Increment load step
			if ((int) CurrentStep > 1)
				LoadFactor += StepIncrement(simulate);
				
			// Get the force vector
			CurrentStep.Forces = LoadFactor * ForceVector;

			if (!simulate)
				return;
			
			// Calculate the initial increment
			var stiffness = SimplifiedStiffness(GlobalStiffness!, FemInput.ConstraintIndex);
			var f0        = SimplifiedForces(CurrentStep.Forces, FemInput.ConstraintIndex);
			OngoingIteration.DisplacementIncrement = stiffness.Solve(f0);

		}

		/// <summary>
		///		Calculate the current stiffness parameter for defining the load increment.
		/// </summary>
		private double CurrentStiffnessParameter()
		{
			var inc = OngoingIteration.DisplacementIncrement;

			var k = (CurrentStep.Forces.ToRowMatrix() * inc)[0] / (inc.ToRowMatrix() * inc)[0];

			if ((int) CurrentStep > 1)
				return k / _k0;
			
			// Set initial
			_k0 = k;
			
			return 1;
		}

		/// <summary>
		///		Clear the iterations lists.
		/// </summary>
		protected virtual void ClearIterations()
		{
			// Clear iteration list
			if ((int) CurrentStep <= 1 || (int) OngoingIteration < 4)
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
					StopMessage = $"Convergence not reached at step {(int) CurrentStep}";
					return Stop;

				default:
					return
						VerifyConvergence(OngoingIteration.ForceConvergence,        ForceTolerance)          ||
						VerifyConvergence(OngoingIteration.DisplacementConvergence, DisplacementTolerance);
			}
		}

		/// <summary>
		///     Save step results after achieving convergence.
		/// </summary>
		private void SaveLoadStepResults()
		{
			CurrentStep.IsCalculated  = true;
			CurrentStep.Convergence   = OngoingIteration.ForceConvergence;
			CurrentStep.Displacements = DisplacementVector!;
			CurrentStep.Stiffness     = GlobalStiffness!;

			if (!_monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = Length.FromMillimeters(DisplacementVector![_monitoredIndex.Value]);

			// Set to step
			CurrentStep.MonitoredDisplacement = new MonitoredDisplacement(disp, LoadFactor);
		}

		/// <summary>
		///     Execute step by step analysis.
		/// </summary>
		/// <param name="simulate">Set true to execute analysis until convergence is not achieved (structural failure).</param>
		private void StepAnalysis(bool simulate)
		{
			// Initiate first step
			CurrentStep.Number = 1;

			// Step-by-step analysis
			do
			{
				// Iterate
				Iterate(simulate);

				// Verify if convergence was not reached
				if (Stop)
					goto CorrectResults;

				// Set step results
				SaveLoadStepResults();

				// Create step
				_steps.Add(CurrentStep.Clone());

				// Increment step
				CurrentStep.Number++;
			} while (simulate || (int) CurrentStep <= NumberOfSteps);

			CorrectResults:
			CorrectResults();
		}

		/// <summary>
		///		Get the step increment.
		/// </summary>
		/// <inheritdoc cref="StepAnalysis"/>
		private double StepIncrement(bool simulate) => 1D / NumberOfSteps;
		
		/// <summary>
		///     Update displacements.
		/// </summary>
		private void UpdateDisplacements()
		{
			// Increment displacements
			var stiffness = SimplifiedStiffness(GlobalStiffness!, FemInput.ConstraintIndex);
			OngoingIteration.DisplacementIncrement = -stiffness.Solve(CurrentSolution.ResidualForces);
			DisplacementVector                     += OngoingIteration.DisplacementIncrement;

			// Update displacements in grips and elements
			FemInput.Grips.SetDisplacements(DisplacementVector);
			FemInput.UpdateDisplacements();
		}

		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="convergence"> Calculated convergence. </param>
		/// <param name="tolerance">The required tolerance.</param>
		/// <seealso cref="ForceTolerance"/>
		/// <seealso cref="DisplacementTolerance"/>
		private bool VerifyConvergence(double convergence, double tolerance) => convergence <= tolerance && (int) OngoingIteration >= MinIterations;

		#endregion

	}
}