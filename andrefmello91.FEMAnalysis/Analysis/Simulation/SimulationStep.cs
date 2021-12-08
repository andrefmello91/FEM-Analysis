using System;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis;

/// <summary>
///     Class for simulation load steps
/// </summary>
public class SimulationStep : LoadStep
{

	#region Fields

	/// <summary>
	///     The initial load factor of this step.
	/// </summary>
	private readonly double _initialLoadFactor;

	#endregion

	#region Properties

	/// <summary>
	///     The Arc-Length calculated for this load step.
	/// </summary>
	public double ArcLength { get; private set; }

	/// <summary>
	///     The desired number of iterations for achieving convergence for this load step.
	/// </summary>
	/// <remarks>
	///     Default : 5
	/// </remarks>
	public int DesiredIterations { get; set; } = 10;

	/// <inheritdoc />
	public override double LoadFactor => ((SimulationIteration) CurrentIteration).LoadFactor;

	/// <inheritdoc />
	public override double LoadFactorIncrement => LoadFactor - _initialLoadFactor;

	/// <summary>
	///     The required number of iterations for achieving convergence for this load step.
	/// </summary>
	public int RequiredIterations => Iterations.Count(i => i.Number > 0);

	#endregion

	#region Constructors

	/// <inheritdoc />
	internal SimulationStep(ForceVector fullForceVector, double loadFactor, AnalysisParameters parameters, int number = 0)
		: base(fullForceVector, loadFactor, parameters, number, true) =>
		_initialLoadFactor = loadFactor;

	/// <inheritdoc />
	internal SimulationStep(int number, ForceVector fullForceVector, double loadFactor, DisplacementVector initialDisplacements, StiffnessMatrix stiffness, AnalysisParameters parameters)
		: base(number, fullForceVector, loadFactor, initialDisplacements, stiffness, parameters, true) =>
		_initialLoadFactor = loadFactor;

	#endregion

	#region Methods

	/// <inheritdoc cref="LoadStep.FromLastStep" />
	public static SimulationStep FromLastStep(SimulationStep lastStep)
	{
		var newStep = (SimulationStep) From(lastStep.FullForceVector, lastStep.LoadFactor, lastStep.FinalDisplacements, lastStep.Stiffness, lastStep.Parameters, lastStep.Number + 1, true);

		// Set desired iterations
		newStep.DesiredIterations = lastStep.DesiredIterations;

		// Update arc length
		newStep.CalculateArcLength(lastStep);

		return newStep;
	}

	/// <summary>
	///     Do the initial load step of a nonlinear simulation procedure.
	/// </summary>
	/// <inheritdoc cref="LoadStep.InitialStep" />
	/// <returns>
	///     The initial <see cref="SimulationStep" />.
	/// </returns>
	internal static SimulationStep InitialStep(IFEMInput femInput, AnalysisParameters parameters)
	{
		var step = (SimulationStep) From(femInput, 0, parameters, 1, true);

		// Set initial increment
		var iteration = (SimulationIteration) step.CurrentIteration;
		iteration.LoadFactorIncrement = 0.3;

		// Get the initial stiffness and force vector simplified
		iteration.Stiffness = femInput.AssembleStiffness();

		// Set initial residual
		var intForces = ForceVector.Zero(femInput.NumberOfDoFs);
		iteration.UpdateForces(step.Forces, intForces);

		// Calculate the initial increments
		var dUr = DisplacementVector.Zero(femInput.NumberOfDoFs);
		var dUf = iteration.Stiffness.Solve(step.FullForceVector);
		iteration.IncrementDisplacements(dUr, dUf, true);

		// Calculate arc length
		step.ArcLength = InitialArcLenght(iteration);

		return step;
	}

	/// <summary>
	///     Get the initial arc lenght.
	/// </summary>
	private static double InitialArcLenght(SimulationIteration initialIteration) =>
		initialIteration.LoadFactorIncrement * (initialIteration.IncrementFromExternal.ToRowMatrix() * (Vector<double>) initialIteration.IncrementFromExternal)[0].Sqrt();

	/// <summary>
	///     Get the accumulated load factor increment from the beginning of this step until the <paramref name="finalIndex" />
	///     iteration.
	/// </summary>
	/// <param name="finalIndex">The required final index to get the increment.</param>
	public double AccumulatedLoadFactorIncrement(Index finalIndex)
	{
		var iterations = Iterations.Where(i => i.Number > 0).ToArray();

		double accL;

		try
		{
			accL = ((SimulationIteration) iterations[finalIndex]).LoadFactor - _initialLoadFactor;
		}
		catch
		{
			accL = 0;
		}

		return accL;
	}

	/// <inheritdoc />
	public override void Iterate(IFEMInput femInput)
	{
		if (Converged)
			return;

		// Initiate first iteration
		foreach (var iteration in Iterations)
			iteration.Number = 0;

		// Iterate
		do
		{
			// Update stiffness
			UpdateStiffness();

			// Add iteration
			NewIteration(true);

			// Update displacements
			UpdateDisplacements();

			// Calculate increment
			IterationIncrement();

			// Update displacements in grips and elements
			((SimulationIteration) CurrentIteration).UpdateDisplacements();

			// Update elements and forces
			UpdateElements(femInput);
			UpdateForces(femInput);

			// Calculate convergence
			CurrentIteration.CalculateConvergence(Forces, FirstIteration.DisplacementIncrement);

		} while (!IterativeStop());
	}

	/// <summary>
	///     Calculate the arc length.
	/// </summary>
	/// <param name="lastStep">The last calculated load step.</param>
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
				var ds1 = (dU.ToRowMatrix() * (Vector<double>) dU)[0].Sqrt();
				ArcLength = ds1 * DesiredIterations / lastStep.RequiredIterations;
				return;
		}
	}

	/// <summary>
	///     Calculate and set the load increment for the current iteration.
	/// </summary>
	private void IterationIncrement()
	{
		var curIt = (SimulationIteration) CurrentIteration;
		var dUf   = curIt.IncrementFromExternal;
		var dUr   = curIt.IncrementFromResidual;
		var dS    = ArcLength;

		// Get accumulated increment until last iteration
		var deltaU = AccumulatedDisplacementIncrement(^2);

		// Calculate coefficients
		var ds2       = dS * dS;
		var a1        = (dUf.ToRowMatrix() * (Vector<double>) dUf)[0];
		var dUrPlusDu = dUr + deltaU;
		var a2        = (dUrPlusDu.ToRowMatrix() * (Vector<double>) dUf)[0];
		var a3        = (dUrPlusDu.ToRowMatrix() * (Vector<double>) dUrPlusDu)[0] - ds2;

		// Calculate roots
		var (r1, r2) = FindRoots.Quadratic(a3, 2 * a2, a1);
		var d1 = r1.Real;
		var d2 = r2.Real;

		if (curIt <= 1)
		{
			// Calculate stiffness determinant
			var det = curIt.Stiffness.Simplified().Determinant() >= 0
				? 1
				: -1;

			// Select the same sign
			curIt.LoadFactorIncrement = d1 / det >= 0
				? d1
				: d2;

			return;
		}

		// Calculate increments and products
		var dF1     = curIt.LoadFactor + d1;
		var dF2     = curIt.LoadFactor + d2;
		var dU1     = -curIt.Stiffness.Solve((ForceVector) (curIt.InternalForces - dF1 * FullForceVector)) + d1 * dUf;
		var dU2     = -curIt.Stiffness.Solve((ForceVector) (curIt.InternalForces - dF2 * FullForceVector)) + d2 * dUf;
		var deltaU1 = deltaU + dU1;
		var deltaU2 = deltaU + dU2;
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

	/// <inheritdoc cref="LoadStep.UpdateDisplacements" />
	private void UpdateDisplacements()
	{
		var curIt = (SimulationIteration) CurrentIteration;
		var lastIt = Iterations.Count > 1
			? (SimulationIteration) LastIteration
			: curIt;

		// Calculate increment from residual
		var dUr = (DisplacementVector) (-curIt.Stiffness.Solve(lastIt.ResidualForces));

		// Calculate increment from external forces
		var dUf = curIt.Stiffness.Solve(FullForceVector);
		curIt.IncrementDisplacements(dUr, dUf);
	}

	#endregion

}