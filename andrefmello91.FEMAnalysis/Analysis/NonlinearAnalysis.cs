using System.Collections.Generic;
using System.Linq;
using System;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Nonlinear analysis class
	/// </summary>
	public class NonlinearAnalysis : Analysis<INonlinearElement>
	{

		#region Fields

		/// <summary>
		///     The force vector of current load step.
		/// </summary>
		/// <inheritdoc cref="Analysis{TFiniteElement}.ForceVector" />
		private Vector<double> _currentForces;

		/// <summary>
		///     The results of the current iteration.
		/// </summary>
		private IterationResult _currentIteration;
		
		/// <summary>
		///     The results of the last iteration.
		/// </summary>
		private IterationResult _lastIteration;

		/// <summary>
		///     Current load step number.
		/// </summary>
		private int _loadStep;

		/// <summary>
		///     The list of load step results.
		/// </summary>
		private readonly List<LoadStepResult> _loadSteps = new();

		/// <summary>
		///     Field to store the DoF index for monitored displacements.
		/// </summary>
		private int? _monitoredIndex;

		#endregion

		#region Properties

		/// <summary>
		///     Get/set the maximum number of iterations.
		/// </summary>
		public int MaxIterations { get; set; }

		/// <summary>
		///     Get/set the minimum number of iterations.
		/// </summary>
		public int MinIterations { get; set; }

		/// <summary>
		///     Get/set the number of load steps to execute.
		/// </summary>
		public int NumLoadSteps { get; set; }

		/// <summary>
		///     The nonlinear equation solver.
		/// </summary>
		public NonLinearSolver Solver { get; set; }

		/// <summary>
		///     Get/set when to stop analysis.
		/// </summary>
		public bool Stop { get; private set; }

		/// <summary>
		///     Get/set the stop message.
		/// </summary>
		public string StopMessage { get; private set; } = string.Empty;

		/// <summary>
		///     Get/set the convergence tolerance.
		/// </summary>
		public double Tolerance { get; set; }

		/// <summary>
		///     Get current convergence.
		/// </summary>
		private double Convergence
		{
			get
			{
				double
					num = 0,
					den = 1;

				var res = _currentIteration.ResidualForces;
				var f   = _currentForces;

				for (var i = 0; i < res.Count; i++)
				{
					num += res[i] * res[i];
					den += f[i] * f[i];
				}

				return
					num / den;
			}
		}

		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		private bool ConvergenceReached => VerifyConvergence(Convergence);

		/// <summary>
		///     The current load step result.
		/// </summary>
		private LoadStepResult CurrentLoadStep => _loadSteps.Last();

		/// <summary>
		///     The last load step result.
		/// </summary>
		private LoadStepResult LastLoadStep => _loadSteps.Count > 1
			? _loadSteps[_loadSteps.Count - 2]
			: CurrentLoadStep;

		/// <summary>
		///     Get current load factor.
		/// </summary>
		private double LoadFactor => (double) _loadStep / NumLoadSteps;

		#endregion

		#region Constructors

		/// <summary>
		///     Nonlinear analysis constructor.
		/// </summary>
		/// <param name="nonlinearInput">The <see cref="IFEMInput{INonlinearElement}" />.</param>
		/// <param name="solver">The <see cref="NonLinearSolver" /> to use.</param>
		/// <param name="numLoadSteps">The number of load steps to perform (default: 50).</param>
		/// <param name="tolerance">The convergence tolerance (default: 1E-6).</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step (default: 10000).</param>
		/// <param name="minIterations">Minimum number of iterations for each load step (default: 2).</param>
		public NonlinearAnalysis(
			IFEMInput<INonlinearElement> nonlinearInput,
			NonLinearSolver solver = NonLinearSolver.Secant,
			int numLoadSteps = 50,
			double tolerance = 1E-6,
			int maxIterations = 10000,
			int minIterations = 2)
			: base(nonlinearInput)
		{
			Solver        = solver;
			NumLoadSteps  = numLoadSteps;
			Tolerance     = tolerance;
			MaxIterations = maxIterations;
			MinIterations = minIterations;
		}

		#endregion

		#region Methods

		/// <summary>
		///     Execute the analysis.
		/// </summary>
		/// <inheritdoc cref="Initiate" />
		/// <inheritdoc cref="StepAnalysis" />
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis{TFiniteElement}.ForceVector" /> (default: 1).</param>
		public void Execute(int? monitoredIndex = null, bool simulate = false, double loadFactor = 1)
		{
			// Get force vector
			ForceVector = FemInput.ForceVector * loadFactor;

			// Get the initial stiffness and force vector simplified
			base.UpdateStiffness();

			// Initiate lists
			Initiate(monitoredIndex);

			// Analysis by load steps
			StepAnalysis(simulate);

			// Set displacements
			DisplacementVector = CurrentLoadStep.Displacements;
			GlobalStiffness    = CurrentLoadStep.Stiffness;

			// Set Reactions
			FemInput.Grips.SetReactions(GetReactions());
		}

		/// <summary>
		///     Generate an <see cref="FEMOutput" /> from analysis results.
		/// </summary>
		/// <returns>
		///     null if no monitored index was provided.
		/// </returns>
		public FEMOutput? GenerateOutput() =>
			!_monitoredIndex.HasValue
				? null
				: new FEMOutput(_loadSteps.Select(ls => ls.MonitoredDisplacement!.Value).ToList());

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix" /> of current iteration.
		/// </summary>
		/// <param name="simplify">Simplify stiffness?</param>
		protected override void UpdateStiffness(bool simplify = true)
		{
			switch (Solver)
			{
				case NonLinearSolver.Secant:
					// Clone current stiffness
					var kCur = _currentIteration.Stiffness.Clone();

					// Calculate the variation of displacements and residual as vectors
					Vector<double>
						dDisp = _currentIteration.Displacements - _lastIteration.Displacements,
						dRes  = _currentIteration.ResidualForces - _lastIteration.ResidualForces;

					// Increment current stiffness
					var dK = ((dRes - _lastIteration.Stiffness * dDisp) / dDisp.Norm(2)).ToColumnMatrix() * dDisp.ToRowMatrix();

					// Set new values
					_currentIteration.Stiffness = _lastIteration.Stiffness + dK;
					_lastIteration.Stiffness    = kCur;
					
					break;

				// For Newton-Raphson
				default:
					// Update stiffness in elements
					FemInput.Elements.UpdateStiffness();

					// Set new values
					_lastIteration.Stiffness    = _currentIteration.Stiffness.Clone();
					_currentIteration.Stiffness = Extensions.AssembleStiffness(FemInput);

					break;
			}

			// Simplify
			if (simplify)
				Simplify(simplify, false);
		}

		/// <summary>
		///     Correct results from last load step after not achieving convergence.
		/// </summary>
		private void CorrectResults()
		{
			// Set displacements from last load step
			FemInput.Grips.SetDisplacements(LastLoadStep.Displacements);
			FemInput.Elements.UpdateDisplacements();

			// Calculate element forces
			FemInput.Elements.CalculateForces();
		}

		/// <summary>
		///     Update displacements.
		/// </summary>
		private void DisplacementUpdate()
		{
			// Set last displacements
			_lastIteration.Displacements = _currentIteration.Displacements.Clone();

			// Increment displacements
			_currentIteration.Displacements += _currentIteration.Stiffness.Solve(-_currentIteration.ResidualForces);

			// Update displacements in grips and elements
			FemInput.Grips.SetDisplacements(_currentIteration.Displacements);
			FemInput.Elements.UpdateDisplacements();
		}

		/// <summary>
		///     Initiate fields.
		/// </summary>
		/// <param name="monitoredIndex">The index of a degree of freedom to monitor, if wanted.</param>
		private void Initiate(int? monitoredIndex)
		{
			_monitoredIndex = monitoredIndex;

			// Calculate initial displacements
			var lf0 = 1D / NumLoadSteps;
			var di  = GlobalStiffness!.Solve(lf0 * ForceVector);

			// Initiate solution values
			_currentIteration = new IterationResult(di, Vector<double>.Build.Dense(FemInput.NumberOfDoFs), GlobalStiffness!);
			_lastIteration    = _currentIteration.Clone();
		}

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		private void Iterate()
		{
			// Initiate first iteration
			_currentIteration.Number = 1;

			while ((int) _currentIteration <= MaxIterations)
			{
				// Calculate element forces
				FemInput.Elements.CalculateForces();

				// Update residual
				ResidualUpdate();

				// Check convergence or if analysis must stop
				if (ConvergenceReached || StopCheck())
					return;

				// Update stiffness and displacements
				DisplacementUpdate();
				UpdateStiffness();

				// Increase iteration count
				_currentIteration.Number++;
			}
		}

		/// <summary>
		///     Calculate residual force <see cref="Vector{T}" />.
		/// </summary>
		private Vector<double> ResidualForces() => Extensions.InternalForces<INonlinearElement>(FemInput) - _currentForces;

		/// <summary>
		///     Update residual force <see cref="Vector{T}" />.
		/// </summary>
		private void ResidualUpdate()
		{
			// Set new values
			_lastIteration.ResidualForces    = _currentIteration.ResidualForces.Clone();
			_currentIteration.ResidualForces = ResidualForces();
		}

		/// <summary>
		///     Save load step results after achieving convergence.
		/// </summary>
		private void SaveLoadStepResults()
		{
			var curLoadStep = new LoadStepResult(_loadStep, _currentForces, _currentIteration.Displacements, _currentIteration.Stiffness);

			if (_monitoredIndex.HasValue)
			{
				// Get displacement
				var disp = Length.FromMillimeters(_currentIteration.Displacements![_monitoredIndex.Value]);

				// Set to load step
				curLoadStep.MonitoredDisplacement = new MonitoredDisplacement(disp, LoadFactor);
			}

			_loadSteps.Add(curLoadStep);
		}

		/// <summary>
		///     Execute step by step analysis.
		/// </summary>
		/// <param name="simulate">Set true to execute analysis until convergence is not achieved (structural failure).</param>
		private void StepAnalysis(bool simulate)
		{
			// Initiate first load step
			_loadStep = 1;

			while (simulate || _loadStep <= NumLoadSteps)
			{
				// Get the force vector
				_currentForces = LoadFactor * ForceVector;

				// Iterate
				Iterate();

				// Verify if convergence was not reached
				if (Stop)
				{
					CorrectResults();
					return;
				}

				// Set load step results
				SaveLoadStepResults();

				// Increment load step
				_loadStep++;
			}
		}

		/// <summary>
		///     Check if analysis must stop.
		/// </summary>
		private bool StopCheck()
		{
			// Check if one stop condition is reached
			Stop = (int) _currentIteration == MaxIterations || _currentIteration.ResidualForces.ContainsNaNOrInfinity() ||
			       _currentIteration.Displacements.ContainsNaNOrInfinity() || _currentIteration.Stiffness.ContainsNaN();

			// Check if maximum number of iterations is reached
			if (Stop)
				StopMessage = $"Convergence not reached at load step {_loadStep}";

			return Stop;
		}

		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="convergence">
		///     Calculated convergence.
		///     <para>See: <see cref="Convergence" />.</para>
		/// </param>
		private bool VerifyConvergence(double convergence) => convergence <= Tolerance && (int) _currentIteration >= MinIterations;

		#endregion

	}
}