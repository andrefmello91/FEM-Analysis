using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Nonlinear analysis class
	/// </summary>
	public class NonlinearAnalysis : Analysis<IFiniteElement>, IEnumerable<LoadStep>
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

		/// <inheritdoc />
		/// <remarks>
		///     The displacements of current step.
		/// </remarks>
		public override Vector<double> DisplacementVector => CurrentStep.FinalDisplacements;

		/// <inheritdoc />
		/// <remarks>
		///     The stiffness of current step.
		/// </remarks>
		public override Matrix<double> GlobalStiffness => CurrentStep.Stiffness;

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
		/// <param name="nonlinearInput">The <see cref="IFEMInput{INonlinearElement}" />.</param>
		public NonlinearAnalysis(IFEMInput<IFiniteElement> nonlinearInput)
			: this(nonlinearInput, AnalysisParameters.Default)
		{
		}

		/// <summary>
		///     Nonlinear analysis constructor.
		/// </summary>
		/// <param name="nonlinearInput">The <see cref="IFEMInput{INonlinearElement}" />.</param>
		/// <param name="parameters">The analysis parameters.</param>
		public NonlinearAnalysis(IFEMInput<IFiniteElement> nonlinearInput, AnalysisParameters parameters)
			: base(nonlinearInput) =>
			Parameters = parameters;

		#endregion

		#region Methods

		/// <summary>
		///     Calculate the secant stiffness increment.
		/// </summary>
		/// <param name="currentStiffness">The stiffness matrix from current iteration.</param>
		/// <param name="currentDisplacements">The displacement vector from current iteration.</param>
		/// <param name="lastDisplacements">The displacement vector from the last iteration.</param>
		/// <param name="currentResidual">The residual force vector from current iteration.</param>
		/// <param name="lastResidual">The residual force vector from last iteration.</param>
		/// <returns>
		///     <see cref="Matrix{T}" />
		/// </returns>
		public static Matrix<double> SecantIncrement(Matrix<double> currentStiffness, Vector<double> currentDisplacements, Vector<double> lastDisplacements, Vector<double> currentResidual, Vector<double> lastResidual)
		{
			// Calculate the variation of displacements and residual as vectors
			Vector<double>
				dU = currentDisplacements - lastDisplacements,
				dR = currentResidual - lastResidual;

			return
				((dR - currentStiffness * dU) / dU.Norm(2)).ToColumnMatrix() * dU.ToRowMatrix();
		}

		/// <summary>
		///     Calculate the tangent stiffness increment.
		/// </summary>
		/// <param name="currentInternalForces">The internal force vector from current iteration.</param>
		/// <param name="lastInternalForces">The internal force vector from last iteration.</param>
		/// <param name="currentDisplacements">The displacement vector from current iteration.</param>
		/// <param name="lastDisplacements">The displacement vector from the last iteration.</param>
		/// <returns>
		///     <see cref="Matrix{T}" />
		/// </returns>
		public static Matrix<double> TangentIncrement(Vector<double> currentInternalForces, Vector<double> lastInternalForces, Vector<double> currentDisplacements, Vector<double> lastDisplacements)
		{
			// Get variations
			var dF = currentInternalForces - lastInternalForces;
			var dU = currentDisplacements - lastDisplacements;

			return
				dF.ToColumnMatrix() * dU.ToRowMatrix();
		}

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

		/// <summary>
		///     Execute the analysis.
		/// </summary>
		/// <inheritdoc cref="StepAnalysis" />
		/// <param name="monitoredIndex">The index of a degree of freedom to monitor, if wanted.</param>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis{TFiniteElement}.ForceVector" /> (default: 1).</param>
		public void Execute(int? monitoredIndex = null, bool simulate = false, double loadFactor = 1)
		{
			_simulate      = simulate;
			MonitoredIndex = monitoredIndex;

			// Get force vector
			ForceVector = FemInput.ForceVector * loadFactor;

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

				// Create step
				NewStep();
			} while (_simulate || CurrentStep <= Parameters.NumberOfSteps);

			CorrectResults:
			CorrectResults();
		}

		/// <summary>
		///     Get the step increment.
		/// </summary>
		/// <param name="numberOfSteps">The number of load steps.</param>
		public static double StepIncrement(int numberOfSteps) => 1D / numberOfSteps;

		/// <summary>
		///     Create a new load step.
		/// </summary>
		///  <param name="incrementLoad">Increment load of the new step?</param>
		protected void NewStep(bool incrementLoad = true) => Steps.Add(LoadStep.FromLastStep(CurrentStep, incrementLoad));

		#region Interface Implementations

		/// <inheritdoc />
		public IEnumerator<LoadStep> GetEnumerator() => Steps.GetEnumerator();

		/// <inheritdoc />
		IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

		#endregion

		#endregion

	}
}