/*using System;
using andrefmello91.Extensions;
using MathNet.Numerics;
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
			CurrentLoadStep.Number = 1;

			// Step-by-step analysis
			do
			{
				// Increment load step
				if (CurrentLoadStep > 1)
					CurrentLoadStep.IncrementLoad(StepIncrement());
				
				// Iterate
				Iterate();

				// Verify if convergence was not reached
				if (Stop)
					goto CorrectResults;

				// Set step results
				CurrentLoadStep.SetResults(MonitoredIndex);

				// Create step
				Steps.Add(CurrentLoadStep.Clone());

				// Increment step
				CurrentLoadStep.Number++;
			} while (true);

			CorrectResults:
			CorrectResults();
			
		}

		/// <inheritdoc />
		protected override void Iterate()
		{
			// Add iteration
			if (CurrentLoadStep > 1)
				CurrentLoadStep.NewIteration(true);

			// Initiate first iteration
			OngoingIteration.Number = 0;

			// Iterate
			do
			{
				// Add iteration
				CurrentLoadStep.NewIteration(true);

				// Increase iteration count
				OngoingIteration.Number++;
				
				// Update stiffness and displacements
				UpdateDisplacements();
				UpdateStiffness();
				
				// Do first iteration steps
				if (OngoingIteration == 1)
					InitialIteration();

				else
					CurrentLoadStep.IncrementLoad(StepIncrement());
				
				// Calculate element forces
				FemInput.CalculateForces();

				// Update internal forces
				InternalForces = FemInput.AssembleInternalForces();

				// Calculate convergence
				OngoingIteration.CalculateForceConvergence(CurrentLoadStep.Forces);
				OngoingIteration.CalculateDisplacementConvergence(FirstIteration.DisplacementIncrement);

			} while (!IterativeStop());
		}

		/// <summary>
		///		Steps to perform at the initial iteration of a load step.
		/// </summary>
		private void InitialIteration()
		{
			var stiffness = SimplifiedStiffness(OngoingIteration.Stiffness, FemInput.ConstraintIndex);

			switch ((int) CurrentLoadStep)
			{
				// First iteration of first load step
				case 1:
					// Set initial increment
					CurrentLoadStep.IncrementLoad(StepIncrement());

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
					CurrentLoadStep.IncrementLoad(StepIncrement());
					OngoingIteration.IncrementDisplacements(rInc, fInc);
					
					return;
			}
		}

		/// <inheritdoc />
		protected override double StepIncrement()
		{
			var dUf = OngoingIteration.IncrementFromExternal;
			var dUr = OngoingIteration.IncrementFromResidual;
			var dS  = OngoingIteration.ArcLength;
			
			switch ((int) OngoingIteration)
			{
				// First iteration of first load step
				case 1 when CurrentLoadStep == 1:
					return base.StepIncrement();
				
				// First iteration of any load step except the first
				case 1:

					// Get the sign
					var i = OngoingIteration.StiffnessParameter >= 0
						? 1
						: -1;

					return
						i * dS * (dUf.ToRowMatrix() * dUf)[0].Pow(-0.5);
				
				// Any other iteration
				default:
					
					// Get accumulated increment until last iteration
					var deltaU = CurrentLoadStep.AccumulatedDisplacementIncrement(^2);
					
					// Calculate coefficients
					var a1        = (dUf.ToRowMatrix() * dUf)[0];
					var dUrPlusDu = dUr + deltaU;
					var a2        = (dUrPlusDu.ToRowMatrix() * dUf)[0];
					var a3        = (dUrPlusDu.ToRowMatrix() * dUrPlusDu)[0] - dS * dS;
					
					// Calculate roots
					var (r1, r2) = FindRoots.Quadratic(a3, 2 * a2, a1);
					var d1    = r1.Real;
					var d2    = r2.Real;
					
					// Choose value
					var deltaU1 = deltaU + dUr + d1 * dUf;
					var deltaU2 = deltaU + dUr + d2 * dUf;
					var p1      = deltaU * deltaU1;
					var p2      = deltaU * deltaU2;
					
					// Check products
					switch (p1)
					{
						case >= 0 when p2 < 0:
							return d1;

						case < 0 when p2 >= 0:
							return d2;

						default:
						{
							// Calculate coefficients
							var dUrPlusDu1 = dUr + deltaU1;
							var dUrPlusDu2 = dUr + deltaU2;
							var a21        = (dUrPlusDu1.ToRowMatrix() * dUf)[0];
							var a22        = (dUrPlusDu2.ToRowMatrix() * dUf)[0];
							var a31        = (dUrPlusDu1.ToRowMatrix() * dUrPlusDu1)[0] - dS * dS;
							var a32        = (dUrPlusDu2.ToRowMatrix() * dUrPlusDu2)[0] - dS * dS;

							return -a31 / a21 <= -a32 / a22
								? d1
								: d2;
						}
					}
			}
		}

		/// <summary>
		///		Calculate the displacement increment.
		/// </summary>
		private void IncrementDisplacements()
		{
			var stiffness = SimplifiedStiffness(OngoingIteration.Stiffness, FemInput.ConstraintIndex);

			switch ((int) CurrentLoadStep)
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
			switch ((int) CurrentLoadStep)
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

			var k = (CurrentLoadStep.Forces.ToRowMatrix() * inc)[0] / (inc.ToRowMatrix() * inc)[0];

			OngoingIteration.StiffnessParameter = k / FirstIteration.StiffnessParameter;
		}
	}
}*/
