using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Nonlinear analysis class
	/// </summary>
	public class NonlinearAnalysis : Analysis, IEnumerable<LoadStep>
	{

		#region Fields

		/// <summary>
		///     The list of step results.
		/// </summary>
		protected readonly List<LoadStep> Steps = new();

		/// <summary>
		///     Set true to execute analysis until convergence is not achieved (structural failure).
		/// </summary>
		private bool _simulate;

		/// <summary>
		///     Field to store the DoF index for monitored displacements.
		/// </summary>
		protected int? MonitoredIndex;

		#endregion

		#region Properties

		/// <inheritdoc cref="List{T}.this[int]" />
		public LoadStep this[int index] => Steps[index];

		/// <inheritdoc cref="List{T}.this[int]" />
		public LoadStep this[Index index] => Steps[index];

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

		/// <summary>
		///     The current step result.
		/// </summary>
		protected LoadStep CurrentStep => Steps[^1];

		/// <summary>
		///     The last step result.
		/// </summary>
		protected LoadStep LastStep => Steps.Count > 1
			? Steps[^2]
			: CurrentStep;

		#endregion

		#region Constructors

		/// <summary>
		///     Nonlinear analysis constructor with default parameters.
		/// </summary>
		/// <param name="nonlinearInput">The finite element input>.</param>
		public NonlinearAnalysis(IFEMInput nonlinearInput)
			: this(nonlinearInput, AnalysisParameters.Default)
		{
		}

		/// <summary>
		///     Nonlinear analysis constructor.
		/// </summary>
		/// <inheritdoc cref="NonlinearAnalysis(IFEMInput)" />
		/// <param name="parameters">The analysis parameters.</param>
		public NonlinearAnalysis(IFEMInput nonlinearInput, AnalysisParameters parameters)
			: base(nonlinearInput) =>
			Parameters = parameters;

		#endregion

		#region Methods

		/// <summary>
		///     Get the step increment.
		/// </summary>
		/// <param name="numberOfSteps">The number of load steps.</param>
		public static double StepIncrement(int numberOfSteps) => 1D / numberOfSteps;

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
		///     Execute the analysis.
		/// </summary>
		/// <inheritdoc cref="StepAnalysis" />
		/// <param name="monitoredIndex">The index of a degree of freedom to monitor, if wanted.</param>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis.Forces" /> (default: 1).</param>
		public void Execute(int? monitoredIndex = null, bool simulate = false, double loadFactor = 1)
		{
			_simulate      = simulate;
			MonitoredIndex = monitoredIndex;

			// Get force vector
			if (!loadFactor.Approx(1))
				Forces = (ForceVector) (Forces * loadFactor);

			// Analysis by steps
			StepAnalysis();

			// Set Reactions
			FemInput.Grips.SetReactions(GetReactions());
		}

		/// <summary>
		///     Generate an <see cref="FEMOutput" /> from analysis results.
		/// </summary>
		/// <returns>
		///     null if no monitored index was provided.
		/// </returns>
		public FEMOutput GenerateOutput() => new(Steps);

		/// <summary>
		///     Correct results from last step after not achieving convergence.
		/// </summary>
		protected void CorrectResults()
		{
			StopMessage = $"Convergence not reached at {CurrentStep}";

			// Set displacements from last (current now) step
			FemInput.Grips.SetDisplacements(LastStep.FinalDisplacements);
			FemInput.UpdateDisplacements();

			// Calculate element forces
			FemInput.CalculateForces();
		}

		/// <summary>
		///     Initiate fields.
		/// </summary>
		protected virtual void InitialStep()
		{
			// Do the initial step
			var step = LoadStep.InitialStep(FemInput, Parameters, _simulate);

			// Initiate lists solution values
			Steps.Clear();
			Steps.Add(step);
		}

		/// <summary>
		///     Create a new load step.
		/// </summary>
		/// <param name="incrementLoad">Increment load of the new step?</param>
		protected void NewStep(bool incrementLoad = true) => Steps.Add(LoadStep.FromLastStep(CurrentStep, incrementLoad));

		/// <summary>
		///     Execute step by step analysis.
		/// </summary>
		protected virtual void StepAnalysis()
		{
			// Initiate solution values
			InitialStep();

			// Step-by-step analysis
			do
			{
				// Iterate
				CurrentStep.Iterate(FemInput);

				// Verify if convergence was not reached
				if (Stop)
					goto CorrectResults;

				// Set step results
				CurrentStep.SetResults(MonitoredIndex);

				// break;

				if (!_simulate && CurrentStep >= Parameters.NumberOfSteps)
					return;

				// Create step
				NewStep();
			} while (_simulate || CurrentStep <= Parameters.NumberOfSteps);

			CorrectResults:
			CorrectResults();
		}

		#region Interface Implementations

		/// <inheritdoc />
		IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

		/// <inheritdoc />
		public IEnumerator<LoadStep> GetEnumerator() => Steps.GetEnumerator();

		#endregion

		#endregion

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
		///     Calculate the secant stiffness increment.
		/// </summary>
		/// <inheritdoc cref="TangentIncrement(andrefmello91.FEMAnalysis.IIteration, andrefmello91.FEMAnalysis.IIteration)" />
		private static StiffnessMatrix SecantIncrement(IIteration currentIteration, IIteration lastIteration)
		{
			// Calculate the variation of displacements and residual as vectors
			var dU = currentIteration.Displacements - lastIteration.Displacements;
			var	dR = currentIteration.ResidualForces - lastIteration.ResidualForces;

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
	}
}