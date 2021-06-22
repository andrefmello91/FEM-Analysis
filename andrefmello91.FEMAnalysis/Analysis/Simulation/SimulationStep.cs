using System;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

using static andrefmello91.FEMAnalysis.Analysis<andrefmello91.FEMAnalysis.IFiniteElement>;
using static andrefmello91.FEMAnalysis.NonlinearAnalysis;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Class for simulation load steps
	/// </summary>
	public class SimulationStep : LoadStep
	{
		/// <summary>
		///     The Arc-Length calculated for this load step.
		/// </summary>
		public double ArcLength { get; set; }

		/// <summary>
		///		The required number of iterations for achieving convergence for this load step.
		/// </summary>
		public int RequiredIterations => Iterations.Count(i => i.Number > 0);

		/// <summary>
		///		The desired number of iterations for achieving convergence for this load step.
		/// </summary>
		/// <remarks>
		///		Default : 5
		///	</remarks>
		public int DesiredIterations { get; set; } = 5;

		/// <inheritdoc />
		internal SimulationStep(Vector<double> fullForceVector, double loadFactor, AnalysisParameters parameters, int number = 0)
			: base(fullForceVector, loadFactor, parameters, number, true)
		{
		}

		/// <inheritdoc />
		internal SimulationStep(int number, Vector<double> fullForceVector, double loadFactor, Vector<double> initialDisplacements, Matrix<double> stiffness, AnalysisParameters parameters)
			: base(number, fullForceVector, loadFactor, initialDisplacements, stiffness, parameters, true)
		{
		}

		/// <inheritdoc />
		public override void IncrementLoad(double loadFactorIncrement)
		{
			((SimulationIteration) OngoingIteration).LoadFactorIncrement = loadFactorIncrement;

			base.IncrementLoad(loadFactorIncrement);
		}

		/// <inheritdoc />
		public override void Iterate(IFEMInput<IFiniteElement> femInput)
		{
			// Initiate first iteration
			foreach (var iteration in Iterations)
				iteration.Number = 0;

			// Iterate
			do
			{
				// Add iteration
				NewIteration(true);

				if (OngoingIteration.Number == 1)
				{
					// Do the initial iteration
					InitialIteration(femInput);
						
					if (Number > 1)
						// Go to new iteration
						continue;
				}
				
				// Increment forces
				else
					IncrementLoad(IterationIncrement());
				
				// Update stiffness and displacements
				UpdateDisplacements(femInput);
				UpdateStiffness(femInput);

				// Calculate element forces
				femInput.CalculateForces();

				// Update internal forces
				var extForces = SimplifiedForces(FullForceVector, femInput.ConstraintIndex);
				var intForces = femInput.AssembleInternalForces();
				OngoingIteration.UpdateForces(extForces, intForces);

				// Calculate convergence
				((SimulationIteration) OngoingIteration).CalculateConvergence(FirstIteration.DisplacementIncrement);
				
			} while (!IterativeStop());
		}
		
		/// <summary>
		///		Steps to perform at the initial iteration of a simulation.
		/// </summary>
		private void InitialIteration(IFEMInput<IFiniteElement> femInput)
		{
			var ongIt     = (SimulationIteration) OngoingIteration;
			
			switch (Number)
			{
				// First iteration of first load step
				case 1:
					// Simplify stiffness
					var stiffness = SimplifiedStiffness(OngoingIteration.Stiffness, femInput.ConstraintIndex);

					// Set initial increment
					IncrementLoad(StepIncrement(Parameters.NumberOfSteps));

					// Set initial residual
					var intForces = stiffness * ongIt.Displacements;
					ongIt.UpdateForces(FullForceVector, intForces);

					// Calculate the initial increments
					var dUr = -stiffness.Solve(ongIt.ResidualForces);
					var dUf = stiffness.Solve(FullForceVector);
					ongIt.IncrementDisplacements(dUr, dUf);
					
					// Check increment sign
					ongIt.CalculateStiffnessParameter(this);

					// Calculate arc length
					CalculateArcLength(0);
					
					break;
				
				// First iteration of any load step except the first
				default:
					// Update stiffness
					UpdateStiffness(femInput);
					stiffness = SimplifiedStiffness(OngoingIteration.Stiffness, femInput.ConstraintIndex);

					// Calculate stiffness parameter for increment sign
					ongIt.CalculateStiffnessParameter(this);
					
					// Calculate increments
					var rInc = -stiffness.Solve(CurrentSolution.ResidualForces);
					var fInc =  stiffness.Solve(SimplifiedForces(FullForceVector, femInput.ConstraintIndex));
					
					// Set increments
					IncrementLoad(IterationIncrement());
					ongIt.IncrementDisplacements(rInc, fInc);
					
					break;
			}
			
			// Update displacements in grips and elements
			femInput.Grips.SetDisplacements(ongIt.Displacements);
			femInput.UpdateDisplacements();
		}

		/// <summary>
		///		Get the load increment for an iteration.
		/// </summary>
		private double IterationIncrement()
		{
			var ongIt = (SimulationIteration) OngoingIteration;
			var dUf   = ongIt.IncrementFromExternal;
			var dUr   = ongIt.IncrementFromResidual;
			var dS    = ArcLength;
			
			switch ((int) ongIt)
			{
				// First iteration of first load step
				case 1 when Number == 1:
					return StepIncrement(Parameters.NumberOfSteps);
				
				// First iteration of any load step except the first
				case 1:

					// Get the sign
					var i = ongIt.StiffnessParameter >= 0
						? 1
						: -1;

					return
						i * dS * (dUf.ToRowMatrix() * dUf)[0].Pow(-0.5);
				
				// Any other iteration
				default:
					
					// Get accumulated increment until last iteration
					var deltaU = AccumulatedDisplacementIncrement(^2);
					
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

		/// <inheritdoc />
		protected override void UpdateDisplacements(IFEMInput<IFiniteElement> femInput)
		{
			var ongIt  = (SimulationIteration) OngoingIteration;
			var curSol = (SimulationIteration) CurrentSolution;

			// Increment displacements
			var stiffness = SimplifiedStiffness(ongIt.Stiffness, femInput.ConstraintIndex);
			
			// Calculate increment from residual
			var dUr = -stiffness.Solve(curSol.ResidualForces);
			
			// Calculate increment from external forces
			var dUf = stiffness.Solve(FullForceVector);
			
			ongIt.IncrementDisplacements(dUr, dUf);
				
			// Update displacements in grips and elements
			femInput.Grips.SetDisplacements(ongIt.Displacements);
			femInput.UpdateDisplacements();
		}
		
		/// <summary>
		///		Calculate the arc length.
		/// </summary>
		/// <param name="requiredIterations">The required iterations for achieving convergence in the last load step.</param>
		public void CalculateArcLength(int requiredIterations)
		{
			switch (Number)
			{
				// First iteration of first load step
				case 1:
					var ongIt = (SimulationIteration) OngoingIteration;
					ArcLength = ongIt.LoadFactorIncrement * (ongIt.DisplacementIncrement.ToRowMatrix() * ongIt.DisplacementIncrement)[0].Sqrt();
					return;
				
				// First iteration of any load step except the first
				default:
					var ds0 = ArcLength;
					ArcLength = ds0 * DesiredIterations / requiredIterations;
					return;
			}
		}

	}
}