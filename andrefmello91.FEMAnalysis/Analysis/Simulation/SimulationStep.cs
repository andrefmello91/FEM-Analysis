using System;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;
using static andrefmello91.FEMAnalysis.Analysis<andrefmello91.FEMAnalysis.IFiniteElement>;
using static andrefmello91.FEMAnalysis.NonlinearAnalysis;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Enumeration for initial load increment sign.
	/// </summary>
	public enum IncrementSign
	{
		/// <summary>
		///		For positive increment.
		/// </summary>
		Positive = 1,
		
		/// <summary>
		///		For negative increment.
		/// </summary>
		Negative = -1
	}
	
	/// <summary>
	///		Class for simulation load steps
	/// </summary>
	public class SimulationStep : LoadStep
	{
		/// <summary>
		///		The sign of the initial load increment.
		/// </summary>
		public IncrementSign Sign { get; private set; }
		
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

		/// <summary>
		///		The force vector increment for the.
		/// </summary>
		public Vector<double> ForceIncrement => ((SimulationIteration) CurrentIteration).LoadFactorIncrement * FullForceVector;

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
			
			// Set initial sign
			step.Sign = IncrementSign.Positive;
			
			// Set initial increment
			var iteration = (SimulationIteration) step.CurrentIteration;
			iteration.LoadFactorIncrement = 0.05;
			
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
			// var dUr = -stiffness.Solve(iteration.ResidualForces);
			var dUr = Vector<double>.Build.Dense(femInput.NumberOfDoFs);
			var dUf = stiffness.Solve(fullforces);
			iteration.IncrementDisplacements(dUr, dUf, true);
			
			// Increment load factor
			// step.IncrementLoad(iteration.LoadFactorIncrement);

			// Calculate arc length
			step.ArcLength = InitialArcLenght((SimulationIteration) step.CurrentIteration);

			// // Update displacements in grips and elements
			// femInput.Grips.SetDisplacements(iteration.Displacements);
			// femInput.UpdateDisplacements();
			//
			// // Calculate element forces
			// femInput.CalculateForces();
			//
			// // Update internal forces
			// iteration.UpdateForces(fullForces, femInput.AssembleInternalForces());

			return step;
		}

		///  <inheritdoc cref="LoadStep.FromLastStep"/>
		public static SimulationStep FromLastStep(SimulationStep lastStep, bool incrementLoad = true)
		{
			var newStep = (SimulationStep) From(lastStep.FullForceVector, lastStep.LoadFactor, lastStep.FinalDisplacements, lastStep.Stiffness, lastStep.Parameters, lastStep.Number + 1, true);
			
			// Set initial sign
			newStep.Sign = GetSign(lastStep);
			
			// Set desired iterations
			newStep.DesiredIterations = lastStep.DesiredIterations;

			// Add last step final iteration
			newStep.Iterations.Clear();
			newStep.Iterations.Add(lastStep.CurrentIteration.Clone());
			
			// Update arc length
			newStep.CalculateArcLength(lastStep);
			
			// Increment load
			if (incrementLoad)
				newStep.IncrementLoad(0.05);
			
			return newStep;
		}

		/// <inheritdoc />
		public override void Iterate(IFEMInput<IFiniteElement> femInput)
		{
			// Initiate first iteration
			foreach (var iteration in Iterations)
				iteration.Number = 0;

			// Do the initial iteration
			if (Number > 1)
				InitialIteration(femInput);
			
			if (Converged)
				return;
			
			// Iterate
			do
			{
				// Increment step load factor
				IncrementLoad(((SimulationIteration) CurrentIteration).LoadFactorIncrement);
				
				// Update forces
				UpdateForces(femInput);

				// Add iteration
				NewIteration(true);
				
				// Update stiffness
				UpdateStiffness(femInput);
				
				// Update displacements
				UpdateDisplacements(femInput);
				
				// Update and Increment forces
				IterationIncrement();

				// Update displacements in grips and elements
				((SimulationIteration) CurrentIteration).UpdateDisplacements();
				femInput.Grips.SetDisplacements(CurrentIteration.Displacements);
				femInput.UpdateDisplacements();
				
				// Calculate convergence
				((SimulationIteration) CurrentIteration).CalculateConvergence(FirstIteration.DisplacementIncrement);
				
			} while (!IterativeStop());
		}

		/// <summary>
		///		Get the initial increment sign for the next load step.
		/// </summary>
		/// <param name="lastLoadStep">The last completed load step.</param>
		private static IncrementSign GetSign(SimulationStep lastLoadStep)
		{
			var df = lastLoadStep.Stiffness.Determinant();
			var di = lastLoadStep.FirstIteration.Stiffness.Determinant();

			return (df / di) switch
			{
				>= 0                                                 => lastLoadStep.Sign,
				< 0 when lastLoadStep.Sign is IncrementSign.Positive => IncrementSign.Negative,
				_                                                    => IncrementSign.Positive
			};
		}

		/// <summary>
		///		Steps to perform at the initial iteration of any load step, except the first.
		/// </summary>
		private void InitialIteration(IFEMInput<IFiniteElement> femInput)
		{
			// Add iteration
			NewIteration(true);

			var curIt     = (SimulationIteration) CurrentIteration;

			// Update stiffness
			UpdateStiffness(femInput);
			var stiffness = SimplifiedStiffness(CurrentIteration.Stiffness, femInput.ConstraintIndex);
			
			// Calculate increments
			var rInc = -stiffness.Solve(CurrentIteration.ResidualForces);
			var fInc =  stiffness.Solve(SimplifiedForces(FullForceVector, femInput.ConstraintIndex));
			
			// Set increments
			curIt.IncrementDisplacements(rInc, fInc);
			
			// Increment load
			IterationIncrement();
			IncrementLoad(curIt.LoadFactorIncrement);
			
			// Update displacements in grips and elements
			curIt.UpdateDisplacements();
			femInput.Grips.SetDisplacements(curIt.Displacements);
			femInput.UpdateDisplacements();
			
			// Check convergence
			Converged = CurrentIteration.CheckConvergence(Parameters);
		}

		/// <inheritdoc />
		protected override void UpdateForces(IFEMInput<IFiniteElement> femInput)
		{
			// Calculate element forces
			femInput.CalculateForces();

			// Update internal forces
			var extForces = SimplifiedForces(Forces, femInput.ConstraintIndex);
			var intForces = femInput.AssembleInternalForces();
			CurrentIteration.UpdateForces(extForces, intForces);
		}

		/// <summary>
		///		Calculate and set the load increment for the current iteration.
		/// </summary>
		private void IterationIncrement()
		{
			var curIt = (SimulationIteration) CurrentIteration;
			var dUf   = curIt.IncrementFromExternal;
			var dUr   = curIt.IncrementFromResidual;
			var dS    = ArcLength;
			
			switch ((int) curIt)
			{
				// First iteration of any load step except the first
				case 1 when Number > 1:
					curIt.LoadFactorIncrement = (int) Sign * dS * (dUf.ToRowMatrix() * dUf)[0].Pow(-0.5);
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
							curIt.LoadFactorIncrement = d1;
							return;

						case < 0 when p2 >= 0:
							curIt.LoadFactorIncrement = d2;
							return;

						default:
						{
							// Calculate coefficients
							var dUrPlusDu1 = dUr + deltaU1;
							var dUrPlusDu2 = dUr + deltaU2;
							var a21        = (dUrPlusDu1.ToRowMatrix() * dUf)[0];
							var a22        = (dUrPlusDu2.ToRowMatrix() * dUf)[0];
							var a31        = (dUrPlusDu1.ToRowMatrix() * dUrPlusDu1)[0] - ds2;
							var a32        = (dUrPlusDu2.ToRowMatrix() * dUrPlusDu2)[0] - ds2;
							
							curIt.LoadFactorIncrement = 
								-a31 / a21 < -a32 / a22
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
			var lastIt = (SimulationIteration) LastIteration;

			// Increment displacements
			var stiffness = SimplifiedStiffness(curIt.Stiffness, femInput.ConstraintIndex);
			
			// Calculate increment from residual
			var dUr = -stiffness.Solve(curIt.ResidualForces);
			
			// Calculate increment from external forces
			var f   = SimplifiedForces(FullForceVector, femInput.ConstraintIndex);
			var dUf = stiffness.Solve(f);
			curIt.IncrementDisplacements(dUr, dUf);
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
					var dU  = lastStep.AccumulatedDisplacementIncrement();
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