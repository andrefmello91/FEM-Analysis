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
		/// 	Set true to execute analysis until convergence is not achieved (structural failure).
		/// </summary>
		private bool _simulate;
		
		/// <summary>
		///     The list of step results.
		/// </summary>
		protected readonly List<StepResult> Steps = new();

		/// <summary>
		///     Field to store the DoF index for monitored displacements.
		/// </summary>
		private int? _monitoredIndex;

		/// <summary>
		///		Get the first iteration of the current step.
		/// </summary>
		protected IterationResult FirstIteration => CurrentStep.Find(i => (int) i == 1)!;
		
		#endregion

		#region Properties

		/// <summary>
		///     The results of the ongoing iteration.
		/// </summary>
		protected IterationResult OngoingIteration => CurrentStep.Count > 0 
			? CurrentStep[^1]
			: LastStep[^1];

		/// <summary>
		///     The results of the current solution (last solved iteration [i - 1]).
		/// </summary>
		protected IterationResult CurrentSolution => CurrentStep.Count > 1
			? CurrentStep[^2]
			: LastStep[^1];

		/// <summary>
		///     The results of the last solution (penultimate solved iteration [i - 2]).
		/// </summary>
		protected IterationResult LastSolution => CurrentStep.Count > 2
			? CurrentStep[^3]
			: LastStep[^2];

		/// <inheritdoc />
		/// <remarks>
		///     The displacements of current step.
		/// </remarks>
		public override Vector<double>? DisplacementVector
		{
			get => CurrentStep.Displacements;
			protected set
			{
				if (value is null)
					return;

				CurrentStep.Displacements = value;
			}
		}

		/// <inheritdoc />
		/// <remarks>
		///     The stiffness of current step.
		/// </remarks>
		public override Matrix<double>? GlobalStiffness
		{
			get => CurrentStep.Stiffness;
			protected set
			{
				if (value is null)
					return;

				CurrentStep.Stiffness = value;
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
		public bool Stop { get; protected set; }

		/// <summary>
		///     Get/set the stop message.
		/// </summary>
		public string StopMessage { get; protected set; } = string.Empty;

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
		protected StepResult CurrentStep => Steps[^1];
		
		/// <summary>
		///     The last step result.
		/// </summary>
		protected StepResult LastStep => Steps.Count > 1
			? Steps[^2]
			: CurrentStep;

		/// <summary>
		///     Get/set the internal force vector of current iteration.
		/// </summary>
		protected virtual Vector<double> InternalForces
		{
			get => OngoingIteration.InternalForces;
			set => OngoingIteration.UpdateForces(SimplifiedForces(CurrentStep.Forces, FemInput.ConstraintIndex), value);
		}

		/// <summary>
		///     Get the residual force vector of current iteration.
		/// </summary>
		protected Vector<double> ResidualForces => OngoingIteration.ResidualForces;

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
		/// <inheritdoc cref="StepAnalysis" />
		/// <param name="monitoredIndex">The index of a degree of freedom to monitor, if wanted.</param>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis{TFiniteElement}.ForceVector" /> (default: 1).</param>
		public void Execute(int? monitoredIndex = null, bool simulate = false, double loadFactor = 1)
		{
			_simulate       = simulate;
			_monitoredIndex = monitoredIndex;

			// Get force vector
			ForceVector = FemInput.ForceVector * loadFactor;
			
			// Analysis by steps
			StepAnalysis();

			// Set Reactions
			FemInput.Grips.SetReactions(GetReactions());
		}

		/// <summary>
		///     Generate an <see cref="FEMOutput" /> from analysis results.
		/// </summary>
		/// <returns>
		///     null if no monitored index was provided.
		/// </returns>
		public FEMOutput GenerateOutput() => new(Steps);

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix{T}" /> of current iteration.
		/// </summary>
		protected override void UpdateStiffness()
		{
			switch (Solver)
			{
				case NonLinearSolver.Secant:
					// Increment current stiffness
					OngoingIteration.Stiffness += SecantIncrement(CurrentSolution.Stiffness, CurrentSolution.Displacements, LastSolution.Displacements, CurrentSolution.ResidualForces, LastSolution.ResidualForces);
					break;

				// For Newton-Raphson
				case NonLinearSolver.NewtonRaphson:
				case NonLinearSolver.ModifiedNewtonRaphson when (int) OngoingIteration == 1:
					// Update stiffness in elements
					FemInput.UpdateStiffness();
					
					// Set new values
					OngoingIteration.Stiffness = FemInput.AssembleStiffness();
					
					break;
				
				default:
					return;
			}
		}

		/// <summary>
		///     Correct results from last step after not achieving convergence.
		/// </summary>
		protected void CorrectResults()
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
		protected virtual void InitialStep()
		{
			// Get initial load factor
			var lf0 = StepIncrement();
			
			// Initiate lists solution values
			Steps.Clear();
			Steps.Add(new StepResult(lf0 * ForceVector, 1) { LoadFactor = lf0 });
			
			CurrentStep.Add(new IterationResult(FemInput.NumberOfDoFs));

			// Get the initial stiffness and force vector simplified
			OngoingIteration.Stiffness = FemInput.AssembleStiffness();
			var stiffness = SimplifiedStiffness(OngoingIteration.Stiffness, FemInput.ConstraintIndex);

			// Calculate initial displacements
			var fi = SimplifiedForces(CurrentStep.Forces, FemInput.ConstraintIndex);
			OngoingIteration.Displacements = stiffness.Solve(fi);

			// Update displacements in grips and elements
			FemInput.Grips.SetDisplacements(OngoingIteration.Displacements);
			FemInput.UpdateDisplacements();

			// Calculate element forces
			FemInput.CalculateForces();

			// Update internal forces
			InternalForces = FemInput.AssembleInternalForces();
		}

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		protected virtual void Iterate()
		{
			// Add iteration
			if (CurrentStep > 1)
				CurrentStep.Add(IterationResult.FromStepResult(LastStep));

			// Initiate first iteration
			OngoingIteration.Number = 0;

			// Iterate
			do
			{
				// Add iteration
				CurrentStep.Add(OngoingIteration.Clone());

				// Increase iteration count
				OngoingIteration.Number++;

				// Update stiffness and displacements
				UpdateDisplacements();
				UpdateStiffness();

				// Calculate element forces
				FemInput.CalculateForces();

				// Update internal forces
				InternalForces = FemInput.AssembleInternalForces();

				// Calculate convergence
				OngoingIteration.CalculateForceConvergence(CurrentStep.Forces);
				OngoingIteration.CalculateDisplacementConvergence(FirstIteration.DisplacementIncrement);

			} while (!IterativeStop());
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
		protected bool IterativeStop()
		{
			// Check if one stop condition is reached
			Stop = OngoingIteration >= MaxIterations                      || ResidualForces.ContainsNaNOrInfinity() ||
			       OngoingIteration.Displacements.ContainsNaNOrInfinity() || OngoingIteration.Stiffness.ContainsNaN();

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
		protected void SaveLoadStepResults()
		{
			CurrentStep.IsCalculated  = true;
			CurrentStep.Convergence   = OngoingIteration.ForceConvergence;
			CurrentStep.Displacements = OngoingIteration.Displacements;
			CurrentStep.Stiffness     = OngoingIteration.Stiffness;

			if (!_monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = Length.FromMillimeters(CurrentStep.Displacements[_monitoredIndex.Value]);

			// Set to step
			CurrentStep.MonitoredDisplacement = new MonitoredDisplacement(disp, CurrentStep.LoadFactor);
		}

		/// <summary>
		///     Execute step by step analysis.
		/// </summary>
		protected virtual void StepAnalysis()
		{
			// Initiate solution values
			InitialStep();

			// Initiate first step
			CurrentStep.Number = 1;

			// Step-by-step analysis
			do
			{
				// Increment load step
				if (CurrentStep > 1)
					CurrentStep.LoadFactor += StepIncrement();
				
				// Get the force vector
				CurrentStep.Forces = CurrentStep.LoadFactor * ForceVector;

				// Iterate
				Iterate();

				// Verify if convergence was not reached
				if (Stop)
					goto CorrectResults;

				// Set step results
				SaveLoadStepResults();

				// Create step
				Steps.Add(CurrentStep.Clone());

				// Increment step
				CurrentStep.Number++;
			} while (_simulate || (int) CurrentStep <= NumberOfSteps);

			CorrectResults:
			CorrectResults();
		}

		/// <summary>
		///		Get the step increment.
		/// </summary>
		/// <inheritdoc cref="StepAnalysis"/>
		protected virtual double StepIncrement() => 1D / NumberOfSteps;
		
		/// <summary>
		///     Update displacements.
		/// </summary>
		protected virtual void UpdateDisplacements()
		{
			// Increment displacements
			var stiffness = SimplifiedStiffness(OngoingIteration.Stiffness, FemInput.ConstraintIndex);
			OngoingIteration.DisplacementIncrement =  -stiffness.Solve(CurrentSolution.ResidualForces);
			OngoingIteration.Displacements        += OngoingIteration.DisplacementIncrement;

			// Update displacements in grips and elements
			FemInput.Grips.SetDisplacements(OngoingIteration.Displacements);
			FemInput.UpdateDisplacements();
		}

		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="convergence"> Calculated convergence. </param>
		/// <param name="tolerance">The required tolerance.</param>
		/// <seealso cref="ForceTolerance"/>
		/// <seealso cref="DisplacementTolerance"/>
		protected bool VerifyConvergence(double convergence, double tolerance) => convergence <= tolerance && OngoingIteration >= MinIterations;

		#endregion

	}
}