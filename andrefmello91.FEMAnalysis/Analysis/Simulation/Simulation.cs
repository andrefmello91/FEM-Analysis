using System;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Simulation class for Arc-Length control analysis.
	/// </summary>
	internal class Simulation : NonlinearAnalysis
	{
		/// <inheritdoc />
		internal Simulation(IFEMInput<IFiniteElement> nonlinearInput, NonLinearSolver solver = NonLinearSolver.NewtonRaphson) : base(nonlinearInput, solver)
		{
		}

		/// <inheritdoc />
		protected override void InitialStep()
		{
			base.InitialStep();
		}
		
		/// <summary>
		///		Get the step increment for a load step.
		/// </summary>
		/// <param name="loadStep">The current load step.</param>
		public static double StepIncrement(LoadStep loadStep)
		{
			var ongIt = (SimulationIteration) loadStep.OngoingIteration;
			var dUf   = ongIt.IncrementFromExternal;
			var dUr   = ongIt.IncrementFromResidual;
			var dS    = ongIt.ArcLength;
			
			switch ((int) ongIt)
			{
				// First iteration of first load step
				case 1 when loadStep == 1:
					return StepIncrement(loadStep.Parameters.NumberOfSteps);
				
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
					var deltaU = loadStep.AccumulatedDisplacementIncrement(^2);
					
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
		public static void CalculateArcLength(LoadStep step)
		{
			if (step.OngoingIteration is not SimulationIteration ongIt)
				return;
			
			switch ((int) step)
			{
				// First iteration of first load step
				case 1 when ongIt <= 1:
					ongIt.ArcLength = ongIt.LoadFactorIncrement * (ongIt.DisplacementIncrement.ToRowMatrix() * ongIt.DisplacementIncrement)[0].Sqrt();
					return;
				
				// First iteration of any load step except the first
				default:
					var ds0 = ((SimulationIteration) step.Last()).ArcLength;
					ongIt.ArcLength = ds0 * step.DesiredIterations / step.RequiredIterations;
					return;
			}
		}
	}
}
