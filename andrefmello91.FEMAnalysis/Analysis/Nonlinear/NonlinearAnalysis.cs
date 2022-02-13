using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

namespace andrefmello91.FEMAnalysis;

/// <summary>
///     Nonlinear analysis class
/// </summary>
public class NonlinearAnalysis : Analysis, IEnumerable<LoadStep>
{

	/// <summary>
	///     The list of step results.
	/// </summary>
	protected readonly List<LoadStep> Steps = new();

	/// <summary>
	///     Set true to execute analysis until convergence is not achieved (structural failure).
	/// </summary>
	private readonly bool _simulate;

	/// <summary>
	///     Field to store the DoF index for monitored displacements.
	/// </summary>
	protected int? MonitoredIndex;

	/// <inheritdoc cref="List{T}.this[int]" />
	public LoadStep this[int index] => Steps[index];

	/// <inheritdoc cref="List{T}.this[int]" />
	public LoadStep this[Index index] => Steps[index];

	/// <summary>
	///     The current step result.
	/// </summary>
	public LoadStep CurrentStep => Steps[^1];

	/// <summary>
	///     The last step result.
	/// </summary>
	public LoadStep LastStep => Steps.Count > 1
		? Steps[^2]
		: CurrentStep;

	/// <summary>
	///     The analysis parameters.
	/// </summary>
	public AnalysisParameters Parameters { get; }

	/// <summary>
	///     Get/set when to stop analysis.
	/// </summary>
	public bool Stop => CurrentStep.Stop;

	/// <summary>
	///     Get/set the stop message.
	/// </summary>
	public string StopMessage { get; protected set; } = string.Empty;

	/// <inheritdoc />
	public override event EventHandler? AnalysisAborted;

	/// <inheritdoc />
	public override event EventHandler? AnalysisComplete;

	/// <summary>
	///     Event to execute when the current load step is aborted.
	/// </summary>
	/// <remarks>
	///     The aborted step is passed to event args.
	/// </remarks>
	public event EventHandler<StepEventArgs>? StepAborted;

	/// <summary>
	///     Event to execute when the current load step converges.
	/// </summary>
	/// <remarks>
	///     The converged step is passed to event args.
	/// </remarks>
	public event EventHandler<StepEventArgs>? StepConverged;

	/// <summary>
	///     Nonlinear analysis constructor with default parameters.
	/// </summary>
	/// <param name="nonlinearInput">The finite element input>.</param>
	/// <param name="monitoredIndex">The index of a degree of freedom to monitor, if wanted.</param>
	/// <param name="simulate">Execute the analysis until convergence is not reached.</param>
	public NonlinearAnalysis(IFEMInput nonlinearInput, int? monitoredIndex = null, bool simulate = false)
		: this(nonlinearInput, AnalysisParameters.Default, monitoredIndex, simulate)
	{
	}

	/// <summary>
	///     Nonlinear analysis constructor.
	/// </summary>
	/// <param name="parameters">The analysis parameters.</param>
	/// <inheritdoc cref="NonlinearAnalysis(IFEMInput, int?, bool)" />
	public NonlinearAnalysis(IFEMInput nonlinearInput, AnalysisParameters parameters, int? monitoredIndex = null, bool simulate = false)
		: base(nonlinearInput)
	{
		Parameters     = parameters;
		MonitoredIndex = monitoredIndex;
		_simulate      = simulate;
	}

	/// <summary>
	///     Get the step increment.
	/// </summary>
	/// <param name="numberOfSteps">The number of load steps.</param>
	public static double StepIncrement(int numberOfSteps) => 1D / numberOfSteps;

	/// <summary>
	///     Calculate the stiffness increment for nonlinear analysis.
	/// </summary>
	/// <param name="solver">The nonlinear solver.</param>
	/// <returns>
	///     The <see cref="andrefmello91.OnPlaneComponents.StiffnessMatrix" /> increment with current unit.
	/// </returns>
	/// <inheritdoc cref="TangentIncrement" />
	public static StiffnessMatrix StiffnessIncrement(IIteration currentIteration, IIteration lastIteration, NonLinearSolver solver) =>
		solver switch
		{
			NonLinearSolver.Secant => SecantIncrement(currentIteration, lastIteration),
			_                      => TangentIncrement(currentIteration, lastIteration)
		};

	/// <summary>
	///     Calculate the convergence.
	/// </summary>
	/// <param name="numerator">
	///     The residual forces of the current iteration or the initial displacement increment of the
	///     current step.
	/// </param>
	/// <param name="denominator">
	///     The applied forces of the current step or the displacement increment of the current
	///     iteration.
	/// </param>
	internal static double CalculateConvergence(IEnumerable<double> numerator, IEnumerable<double> denominator)
	{
		double
			num = numerator.Sum(n => n * n),
			den = 1 + denominator.Sum(n => n * n);

		return
			num / den;
	}

	/// <inheritdoc cref="CalculateConvergence" />
	internal static double CalculateConvergence<TQuantity, TUnit>(QuantityVector<TQuantity, TUnit> numerator, QuantityVector<TQuantity, TUnit> denominator, IEnumerable<int>? constraintIndexes = null)
		where TQuantity : struct, IQuantity<TUnit>
		where TUnit : Enum
	{
		var unit = numerator.Unit;

		var dd = denominator.Unit.Equals(unit)
			? denominator
			: denominator.Convert(unit);

		return
			CalculateConvergence(numerator.Values, dd.Values);
	}

	/// <summary>
	///     Calculate the secant stiffness increment.
	/// </summary>
	/// <inheritdoc cref="TangentIncrement(andrefmello91.FEMAnalysis.IIteration, andrefmello91.FEMAnalysis.IIteration)" />
	private static StiffnessMatrix SecantIncrement(IIteration currentIteration, IIteration lastIteration)
	{
		// Calculate the variation of displacements and residual as vectors
		var dU = currentIteration.Displacements - lastIteration.Displacements;
		var dR = currentIteration.ResidualForces - lastIteration.ResidualForces;

		var unit = dR.Unit.Per(dU.Unit);

		Matrix<double> stiffness = lastIteration.Stiffness.Unit == unit
			? lastIteration.Stiffness
			: lastIteration.Stiffness.Convert(unit);

		var inc = ((dR - stiffness * (Vector<double>) dU) / dU.Norm(2)).ToColumnMatrix() * dU.ToRowMatrix();

		return
			new StiffnessMatrix(inc, unit);
	}

	/// <summary>
	///     Calculate the tangent stiffness increment.
	/// </summary>
	/// <param name="currentIteration">The current iteration.</param>
	/// <param name="lastIteration">The last solved iteration.</param>
	/// <returns>
	///     The <see cref="StiffnessMatrix" /> increment with current unit.
	/// </returns>
	private static StiffnessMatrix TangentIncrement(IIteration currentIteration, IIteration lastIteration)
	{
		// Get variations
		var dF = (ForceVector) (currentIteration.InternalForces - lastIteration.InternalForces);
		var dU = (DisplacementVector) (currentIteration.Displacements - lastIteration.Displacements);

		return
			dF / dU;
	}

	/// <summary>
	///     Execute the full analysis.
	/// </summary>
	/// <param name="loadFactor">The load factor to multiply <see cref="Analysis.Forces" /> (default: 1).</param>
	public void Execute(double loadFactor = 1)
	{
		// Get force vector
		if (!loadFactor.Approx(1))
			Forces = (ForceVector) (Forces * loadFactor);

		// Analysis by steps
		while (true)
		{
			// Add new step
			ExecuteStep();

			if (Stop || !_simulate && CurrentStep >= Parameters.NumberOfSteps)
				return;
		}
	}

	/// <summary>
	///     Add a load step and execute it.
	/// </summary>
	/// <remarks>
	///     Load step is not executed when analysis is aborted or complete.
	/// </remarks>
	public void ExecuteStep()
	{
		// Check conditions
		if (Steps.Any() && (Stop || !_simulate && CurrentStep >= Parameters.NumberOfSteps))
			return;

		// Add new step
		NewStep();

		// Iterate
		CurrentStep.Iterate(FemInput);

		// Verify if convergence was not reached
		if (Stop)
		{
			Invoke(StepAborted, new StepEventArgs(CurrentStep));
			Invoke(AnalysisAborted);
			CorrectResults();
			return;
		}

		// Set step results
		SetStepResults(MonitoredIndex);

		// Check if analysis is done
		if (_simulate || CurrentStep <= Parameters.NumberOfSteps)
			return;

		// Set Reactions
		FemInput.Grips.SetReactions(GetReactions());

		// Invoke analysis complete event
		Invoke(AnalysisComplete);
	}

	/// <summary>
	///     Generate an <see cref="FEMOutput" /> from analysis results.
	/// </summary>
	public FEMOutput GenerateOutput() => new(Steps, FemInput.MonitoredElements.Select(e => e.Monitor!));

	/// <summary>
	///     Correct results from last step after not achieving convergence.
	/// </summary>
	protected virtual void CorrectResults()
	{
		StopMessage = $"Convergence not reached at {CurrentStep}";

		// Set displacements from last (current now) step
		FemInput.Grips.SetDisplacements(LastStep.FinalDisplacements);
		FemInput.UpdateDisplacements();

		// Calculate element forces
		FemInput.CalculateForces();

		// Set Reactions
		FemInput.Grips.SetReactions(GetReactions());
	}

	/// <inheritdoc cref="LoadStep.SetResults" />
	protected virtual void SetStepResults(int? monitoredIndex)
	{
		CurrentStep.SetResults(FemInput, monitoredIndex);
		Invoke(StepConverged, new StepEventArgs(CurrentStep));
	}

	/// <summary>
	///     Create a new load step.
	/// </summary>
	/// <param name="incrementLoad">Increment load of the new step?</param>
	private void NewStep(bool incrementLoad = true) =>
		Steps.Add(Steps.Any()
			? LoadStep.FromLastStep(CurrentStep, incrementLoad)
			: LoadStep.InitialStep(FemInput, Parameters));

	/// <inheritdoc />
	public IEnumerator<LoadStep> GetEnumerator() => Steps.GetEnumerator();

	/// <inheritdoc />
	IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}