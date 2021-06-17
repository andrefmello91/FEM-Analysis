using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis.Simulation
{
	/// <summary>
	///		Simulation class for Arc-Length control analysis.
	/// </summary>
	internal class Simulation : NonlinearAnalysis
	{
		/// <inheritdoc cref="NonlinearAnalysis.OngoingIteration"/>
		protected new SimulationIterationResult OngoingIteration => (SimulationIterationResult) base.OngoingIteration;

		/// <inheritdoc cref="NonlinearAnalysis.CurrentSolution"/>
		protected new SimulationIterationResult CurrentSolution => (SimulationIterationResult) base.CurrentSolution;

		/// <inheritdoc cref="NonlinearAnalysis.LastSolution"/>
		protected new SimulationIterationResult LastSolution => (SimulationIterationResult) base.LastSolution;
		
		/// <inheritdoc cref="NonlinearAnalysis.FirstIteration"/>
		protected new SimulationIterationResult FirstIteration => (SimulationIterationResult) base.FirstIteration;

		/// <inheritdoc />
		protected override Vector<double> InternalForces
		{
			get => OngoingIteration.InternalForces;
			set => OngoingIteration.UpdateForces(SimplifiedForces(ForceVector!, FemInput.ConstraintIndex), value);
		}
		
		/// <inheritdoc />
		internal Simulation(IFEMInput<IFiniteElement> nonlinearInput, NonLinearSolver solver = NonLinearSolver.NewtonRaphson) : base(nonlinearInput, solver)
		{
		}

		/// <inheritdoc />
		protected override void InitialStep()
		{
			base.InitialStep();
		}

		/// <inheritdoc />
		protected override void StepAnalysis()
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
					CurrentStep.IncrementLoad(StepIncrement());
				
				// Iterate
				Iterate();

				// Verify if convergence was not reached
				if (Stop)
					goto CorrectResults;

				// Set step results
				CurrentStep.SetResults(MonitoredIndex);

				// Create step
				Steps.Add(CurrentStep.Clone());

				// Increment step
				CurrentStep.Number++;
			} while (true);

			CorrectResults:
			CorrectResults();
			
		}

		/// <inheritdoc />
		protected override void Iterate()
		{
			// Add iteration
			if (CurrentStep > 1)
				CurrentStep.NewIteration(true);

			// Initiate first iteration
			OngoingIteration.Number = 0;

			// Iterate
			do
			{
				// Add iteration
				CurrentStep.NewIteration(true);

				// Increase iteration count
				OngoingIteration.Number++;
				
				// Update stiffness and displacements
				UpdateDisplacements();
				UpdateStiffness();
				
				// Do first iteration steps
				if (OngoingIteration == 1)
					InitialIteration();

				else
					CurrentStep.IncrementLoad(StepIncrement());
				
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
		///		Steps to perform at the initial iteration of a load step.
		/// </summary>
		private void InitialIteration()
		{
			var stiffness = SimplifiedStiffness(OngoingIteration.Stiffness, FemInput.ConstraintIndex);

			switch ((int) CurrentStep)
			{
				// First iteration of first load step
				case 1:
					// Set initial increment
					CurrentStep.IncrementLoad(base.StepIncrement());

					// Set initial residual
					var intForces = stiffness * OngoingIteration.Displacements;
					OngoingIteration.UpdateForces(ForceVector!, intForces);

					// Calculate the initial increment
					IncrementDisplacements();

					// Calculate arc length
					CalculateArcLength();
					
					return;
				
				// First iteration of any load step except the first
				default:
					// Check increment sign
					CalculateStiffnessParameter();
					
					// Calculate increments
					var rInc = -stiffness.Solve(CurrentSolution.ResidualForces);
					var fInc =  stiffness.Solve(SimplifiedForces(ForceVector!, FemInput.ConstraintIndex));
					
					// Set increments
					CurrentStep.IncrementLoad(StepIncrement());
					OngoingIteration.IncrementDisplacements(rInc, fInc);
					
					return;
			}
		}

		/// <inheritdoc />
		protected override double StepIncrement()
		{
			// Get the sign
			var i = OngoingIteration.StiffnessParameter >= 0
				?  1
				: -1;

			var dUf = OngoingIteration.DisplacementIncrement;
			var dS  = OngoingIteration.ArcLength;

			return
				i * dS * (dUf.ToRowMatrix() * dUf)[0].Pow(-0.5);
		}

		/// <summary>
		///		Calculate the displacement increment.
		/// </summary>
		private void IncrementDisplacements()
		{
			var stiffness = SimplifiedStiffness(OngoingIteration.Stiffness, FemInput.ConstraintIndex);

			switch ((int) CurrentStep)
			{
				// First iteration of first load step
				case 1 when OngoingIteration <= 1:
					// Calculate increment from residual and from full force vector
					var rInc = -stiffness.Solve(OngoingIteration.ResidualForces);
					var fInc = stiffness.Solve(SimplifiedForces(ForceVector!, FemInput.ConstraintIndex));

					// Set increments
					OngoingIteration.IncrementDisplacements(rInc, fInc);
					
					return;
			}
			
		}

		/// <summary>
		///		Calculate the arc length.
		/// </summary>
		private void CalculateArcLength()
		{
			switch ((int) CurrentStep)
			{
				// First iteration of first load step
				case 1 when OngoingIteration == 0:
					OngoingIteration.ArcLength = OngoingIteration.LoadFactorIncrement * (OngoingIteration.DisplacementIncrement.ToRowMatrix() * OngoingIteration.DisplacementIncrement)[0].Sqrt();
					return;
				
				// First iteration of any load step except the first
				default:
					
			}
		}
		
		/// <summary>
		///		Calculate the current stiffness parameter for defining the sign of the load factor increment.
		/// </summary>
		private void CalculateStiffnessParameter()
		{
			if (OngoingIteration <= 1)
			{
				OngoingIteration.StiffnessParameter = 1;
				return;
			}

			var inc = OngoingIteration.DisplacementIncrement;

			var k = (CurrentStep.Forces.ToRowMatrix() * inc)[0] / (inc.ToRowMatrix() * inc)[0];

			OngoingIteration.StiffnessParameter = k / FirstIteration.StiffnessParameter;
		}
	}
}