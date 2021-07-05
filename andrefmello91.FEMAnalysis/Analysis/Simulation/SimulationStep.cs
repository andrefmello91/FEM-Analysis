using System;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

using static andrefmello91.FEMAnalysis.Analysis<andrefmello91.FEMAnalysis.IFiniteElement>;
using static andrefmello91.FEMAnalysis.StiffnessMatrix;
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
		public double ArcLength { get; private set; }

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
		public override double LoadFactorIncrement => AccumulatedLoadFactorIncrement(^1);

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

		/// <summary>
		///     Do the initial load step of a nonlinear simulation procedure.
		/// </summary>
		/// <inheritdoc cref="LoadStep.InitialStep"/>
		/// <returns>
		///     The initial <see cref="SimulationStep" />.
		/// </returns>
		internal static SimulationStep InitialStep(IFEMInput<IFiniteElement> femInput, AnalysisParameters parameters)
		{
			var step = (SimulationStep) From(femInput, StepIncrement(parameters.NumberOfSteps), parameters, 1, true);
			
			// Set initial increment
			var iteration = (SimulationIteration) step.CurrentIteration;
			iteration.LoadFactorIncrement = 0.1;
			
			// Get the initial stiffness and force vector simplified
			iteration.Stiffness = femInput.AssembleStiffness();
			var stiffness = SimplifiedStiffness(iteration.Stiffness, femInput.ConstraintIndex);

			// Calculate initial displacements
			var extForces  = SimplifiedForces(step.Forces, femInput.ConstraintIndex);
			var fullforces = SimplifiedForces(step.FullForceVector, femInput.ConstraintIndex);
			
			// Set initial residual
			var intForces = Vector<double>.Build.Dense(femInput.NumberOfDoFs);
			iteration.UpdateForces(extForces, intForces);

			// Calculate the initial increments
			var dUr = Vector<double>.Build.Dense(femInput.NumberOfDoFs);
			var dUf = stiffness.Solve(fullforces);
			iteration.IncrementDisplacements(dUr, dUf, true);
			
			// Calculate arc length
			step.ArcLength = InitialArcLenght((SimulationIteration) step.CurrentIteration);

			return step;
		}

		///  <inheritdoc cref="LoadStep.FromLastStep"/>
		public static SimulationStep FromLastStep(SimulationStep lastStep)
		{
			var newStep = (SimulationStep) From(lastStep.FullForceVector, lastStep.LoadFactor, lastStep.FinalDisplacements, lastStep.Stiffness, lastStep.Parameters, lastStep.Number + 1, true);
			
			// Set desired iterations
			newStep.DesiredIterations = lastStep.DesiredIterations;

			// Add last step final iteration
			// newStep.Iterations.Clear();
			// newStep.Iterations.Add(lastStep.CurrentIteration.Clone());
			// newStep.CurrentIteration.Number = 0;
			
			// Update arc length
			newStep.CalculateArcLength(lastStep);
			
			return newStep;
		}

		/// <inheritdoc />
		public override void Iterate(IFEMInput<IFiniteElement> femInput)
		{
			if (Converged)
				return;
			
			// Initiate first iteration
			foreach (var iteration in Iterations)
				iteration.Number = 0;
			
			// Iterate
			while (true)
			{
				// Increment step load factor
				IncrementLoad(((SimulationIteration) CurrentIteration).LoadFactorIncrement);
				
				// Update forces
				UpdateForces(femInput);
				
				// Add iteration
				NewIteration(true);
				
				// Update displacements
				UpdateDisplacements(femInput);
				
				// Update and Increment forces
				IterationIncrement(femInput.ConstraintIndex);

				// Update displacements in grips and elements
				((SimulationIteration) CurrentIteration).UpdateDisplacements();
				
				// Calculate convergence
				CurrentIteration.CalculateConvergence(Forces, FirstIteration.DisplacementIncrement);

				// Check convergence or stop criteria
				var stop = IterativeStop();
				
				// Update stiffness
				UpdateStiffness();

				if (stop)
					return;

				// Update elements
				UpdateElements(femInput);
			}
		}

		/// <inheritdoc />
		protected override void UpdateForces(IFEMInput<IFiniteElement> femInput)
		{
			// Update internal forces
			var extForces = SimplifiedForces(Forces, femInput.ConstraintIndex);
			var intForces = femInput.AssembleInternalForces();
			CurrentIteration.UpdateForces(extForces, intForces);
		}

		/// <summary>
		///		Calculate and set the load increment for the current iteration.
		/// </summary>
		private void IterationIncrement(IEnumerable<int> constraintIndex)
		{
			var curIt = (SimulationIteration) CurrentIteration;
			var dUf   = curIt.IncrementFromExternal;
			var dUr   = curIt.IncrementFromResidual;
			var dS    = ArcLength;
			
			switch ((int) curIt)
			{
				// First iteration of any load step except the first
				case 1 when Number > 1:
					
					// Get stiffness determinant sign
					var stiffness = CurrentIteration.Stiffness;
					var sign = stiffness.Determinant() >= 0
						?  1
						: -1;
					
					curIt.LoadFactorIncrement = sign * dS * (dUf.ToRowMatrix() * dUf)[0].Pow(-0.5);

					return;
				
				// Any other iteration
				default:
					// Get accumulated increment until last iteration
					var deltaU = AccumulatedDisplacementIncrement(^2);

					// Calculate coefficients
					var ds2       = dS * dS;
					var a1        = (dUf.ToRowMatrix() * dUf)[0];
					var dUrPlusDu = dUr + deltaU;
					var a2        = (dUrPlusDu.ToRowMatrix() * dUf)[0];
					var a3        = (dUrPlusDu.ToRowMatrix() * dUrPlusDu)[0] - ds2;
					
					// Calculate roots
					var (r1, r2) = FindRoots.Quadratic(a3, 2 * a2, a1);
					var d1      = r1.Real;
					var d2      = r2.Real;
					
					// Calculate increments and products
					var deltaU1 = deltaU + dUr + d1 * dUf;
					var deltaU2 = deltaU + dUr + d2 * dUf;
					var p1      = deltaU * deltaU1;
					var p2      = deltaU * deltaU2;
					
					// Check products
					switch (p1)
					{
						case >= 0 when p2 < 0:
							curIt.LoadFactorIncrement = d1;
							return;
					
						case < 0 when p2 >= 0:
							curIt.LoadFactorIncrement = d2;
							return;
					
						default:
						{
							// Check linear solution
							var ls = -a3 / a2;
							
							// Set the closest value
							curIt.LoadFactorIncrement = (ls - d1).Abs() <= (ls - d2).Abs()
								? d1
								: d2;
							
							return;
						}
					}
			}
		}

		/// <inheritdoc />
		protected override void UpdateDisplacements(IFEMInput<IFiniteElement> femInput)
		{
			var curIt  = (SimulationIteration) CurrentIteration;
			var lastIt = Iterations.Count > 1 
				? (SimulationIteration) LastIteration
				: curIt;

			// Increment displacements
			var stiffness = SimplifiedStiffness(curIt.Stiffness, femInput.ConstraintIndex);
			
			// Calculate increment from residual
			var dUr = -stiffness.Solve(lastIt.ResidualForces);
			
			// Calculate increment from external forces
			var f   = SimplifiedForces(FullForceVector, femInput.ConstraintIndex);
			var dUf = stiffness.Solve(f);
			curIt.IncrementDisplacements(dUr, dUf);
		}

		/// <summary>
		///		Get the accumulated load factor increment from the beginning of this step until the <paramref name="finalIndex"/> iteration.
		/// </summary>
		/// <param name="finalIndex">The required final index to get the increment.</param>
		public double AccumulatedLoadFactorIncrement(Index finalIndex)
		{
			var iterations = Iterations.Where(i => i.Number > 0).ToList();

			return iterations.Count < finalIndex.Value
				? 0
				: iterations.GetRange(0, finalIndex.Value + 1)
					.Cast<SimulationIteration>()
					.Select(i => i.LoadFactorIncrement)
					.Sum();
		}
		
		///  <summary>
		/// 		Calculate the arc length.
		///  </summary>
		///  <param name="lastStep">The last calculated load step.</param>
		private void CalculateArcLength(SimulationStep lastStep)
		{
			switch (Number)
			{
				// First iteration of first load step
				case 1:
					ArcLength = InitialArcLenght((SimulationIteration) CurrentIteration);
					return;
				
				// First iteration of any load step except the first
				default:
					var dU  = lastStep.DisplacementIncrement;
					var ds1 = (dU.ToRowMatrix() * dU)[0].Sqrt();
					ArcLength = ds1 * DesiredIterations / lastStep.RequiredIterations;
					return;
			}
		}

		/// <summary>
		///		Get the initial arc lenght.
		/// </summary>
		private static double InitialArcLenght(SimulationIteration initialIteration) =>
			initialIteration.LoadFactorIncrement * (initialIteration.IncrementFromExternal.ToRowMatrix() * initialIteration.IncrementFromExternal)[0].Sqrt();
	}
}