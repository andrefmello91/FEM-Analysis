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
			} while (true);

			CorrectResults:
			CorrectResults();
			
		}

		/// <inheritdoc />
		protected override void Iterate()
		{
			// Add iteration
			CurrentStep.Add(OngoingIteration.Clone());

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
				
				// Do first iteration steps
				if (OngoingIteration == 1)
					InitialIteration();

				// Update forces
				CurrentStep.Forces = CurrentStep.LoadFactor * ForceVector!;
				
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
			switch ((int) CurrentStep)
			{
				// First iteration of first load step
				case 1:
					// Set initial increment
					OngoingIteration.LoadFactorIncrement = base.StepIncrement();

					// Set initial residual
					var intForces = SimplifiedStiffness(GlobalStiffness!, FemInput.ConstraintIndex) * DisplacementVector!;
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
					
					// Calculate increment
					OngoingIteration.LoadFactorIncrement   = StepIncrement();
					OngoingIteration.IncrementFromResidual = -SimplifiedStiffness(GlobalStiffness!, FemInput.ConstraintIndex).Solve(CurrentSolution.ResidualForces);
					
					// Increment displacements and load factor
					DisplacementVector     += OngoingIteration.DisplacementIncrement;
					CurrentStep.LoadFactor += OngoingIteration.LoadFactorIncrement;
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
			var stiffness = SimplifiedStiffness(GlobalStiffness!, FemInput.ConstraintIndex);

			switch ((int) CurrentStep)
			{
				// First iteration of first load step
				case 1 when OngoingIteration <= 1:
					// Calculate increment from residual and from full force vector
					OngoingIteration.DisplacementIncrement         =  stiffness.Solve(SimplifiedForces(ForceVector!, FemInput.ConstraintIndex));
					OngoingIteration.IncrementFromResidual =  -stiffness.Solve(OngoingIteration.ResidualForces);
					DisplacementVector                             += OngoingIteration.IncrementFromResidual + OngoingIteration.LoadFactorIncrement * OngoingIteration.DisplacementIncrement;
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