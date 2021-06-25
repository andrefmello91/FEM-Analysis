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
			var step = (SimulationStep) From(femInput, 0, parameters, 1, true);
			
			// Set initial sign
			step.Sign = IncrementSign.Positive;
			
			step.IncrementLoad(StepIncrement(parameters.NumberOfSteps));
			
			var iteration = (SimulationIteration) step.CurrentIteration;

			// Get the initial stiffness and force vector simplified
			iteration.Stiffness = femInput.AssembleStiffness();
			var stiffness = SimplifiedStiffness(iteration.Stiffness, femInput.ConstraintIndex);

			// Calculate initial displacements
			var fullForces = SimplifiedForces(step.FullForceVector, femInput.ConstraintIndex);
			
			// Set initial residual
			var intForces = Vector<double>.Build.Dense(femInput.NumberOfDoFs);
			var extForces = SimplifiedForces(step.Forces, femInput.ConstraintIndex);
			iteration.UpdateForces(extForces, intForces);

			// Calculate the initial increments
			// var dUr = -stiffness.Solve(iteration.ResidualForces);
			var dUr = Vector<double>.Build.Dense(femInput.NumberOfDoFs);
			var dUf = stiffness.Solve(fullForces);
			iteration.IncrementDisplacements(dUr, dUf);
					
			// Calculate arc length
			step.CalculateArcLength(0, 0);

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

		///  <summary>
		/// 		Create a load step from the last load step.
		///  </summary>
		///  <param name="lastStep">The last load step.</param>
		public static SimulationStep FromLastStep(SimulationStep lastStep)
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
			newStep.CalculateArcLength(lastStep.ArcLength, lastStep.RequiredIterations);
			
			return newStep;
		}

		/// <inheritdoc />
		public override void IncrementLoad(double loadFactorIncrement)
		{
			((SimulationIteration) CurrentIteration).LoadFactorIncrement = loadFactorIncrement;

			base.IncrementLoad(loadFactorIncrement);
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
			
			// Iterate
			do
			{
				// Add iteration
				NewIteration(true);
				
				// Update stiffness
				UpdateStiffness(femInput);
				
				// Update displacements
				UpdateDisplacements(femInput);
				
				// Update and Increment forces
				UpdateForces(femInput);
				IncrementLoad(IterationIncrement());

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
			var stiffness = SimplifiedStiffness(CurrentIteration.Stiffness, femInput.ConstraintIndex);

			// Update stiffness
			UpdateStiffness(femInput);
			stiffness = SimplifiedStiffness(CurrentIteration.Stiffness, femInput.ConstraintIndex);
			
			// Increment load
			IncrementLoad(IterationIncrement());

			// Calculate increments
			var rInc = -stiffness.Solve(CurrentIteration.ResidualForces);
			var fInc =  stiffness.Solve(SimplifiedForces(FullForceVector, femInput.ConstraintIndex));
			
			// Set increments
			curIt.IncrementDisplacements(rInc, fInc);
			
			// Update displacements in grips and elements
			femInput.Grips.SetDisplacements(curIt.Displacements);
			femInput.UpdateDisplacements();
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
		///		Get the load increment for an iteration.
		/// </summary>
		private double IterationIncrement()
		{
			var curIt = (SimulationIteration) CurrentIteration;
			var dUf   = curIt.IncrementFromExternal;
			var dUr   = curIt.IncrementFromResidual;
			var dS    = ArcLength;
			
			switch ((int) curIt)
			{
				// First iteration of any load step except the first
				case 1 when Number > 1:
					return
						(int) Sign * dS * (dUf.ToRowMatrix() * dUf)[0].Pow(-0.5);
				
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
							var a31        = (dUrPlusDu1.ToRowMatrix() * dUrPlusDu1)[0] - ds2;
							var a32        = (dUrPlusDu2.ToRowMatrix() * dUrPlusDu2)[0] - ds2;
							
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
			var curIt  = (SimulationIteration) CurrentIteration;
			var lastIt = (SimulationIteration) LastIteration;

			// Increment displacements
			var stiffness = SimplifiedStiffness(curIt.Stiffness, femInput.ConstraintIndex);
			
			// Calculate increment from residual
			var dUr = -stiffness.Solve(lastIt.ResidualForces);
			
			// Calculate increment from external forces
			var f   = SimplifiedForces(FullForceVector, femInput.ConstraintIndex);
			var dUf = stiffness.Solve(f);
			
			curIt.IncrementDisplacements(dUr, dUf);
				
			// Update displacements in grips and elements
			femInput.Grips.SetDisplacements(curIt.Displacements);
			femInput.UpdateDisplacements();
		}

		///  <summary>
		/// 		Calculate the arc length.
		///  </summary>
		///  <param name="lastArcLenght">The arc lenght of the last load step.</param>
		///  <param name="requiredIterations">The required iterations for achieving convergence in the last load step.</param>
		private void CalculateArcLength(double lastArcLenght, int requiredIterations)
		{
			switch (Number)
			{
				// First iteration of first load step
				case 1:
					var curIt = (SimulationIteration) CurrentIteration;
					ArcLength = curIt.LoadFactorIncrement * (curIt.IncrementFromExternal.ToRowMatrix() * curIt.IncrementFromExternal)[0].Sqrt();
					return;
				
				// First iteration of any load step except the first
				default:
					ArcLength = lastArcLenght * DesiredIterations / requiredIterations;
					return;
			}
		}

	}
}