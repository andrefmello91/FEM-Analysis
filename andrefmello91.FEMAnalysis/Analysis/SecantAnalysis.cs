using System.Collections.Generic;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Nonlinear analysis class.
	/// </summary>
	public class SecantAnalysis : Analysis
	{

		#region Fields

		/// <summary>
		///     The displacement <see cref="Vector" /> of current iteration
		/// </summary>
		/// <inheritdoc cref="Analysis.DisplacementVector" />
		private Vector<double>? _currentDisplacements;

		/// <summary>
		///     Field to store each iteration force <see cref="Vector" />
		/// </summary>
		/// <inheritdoc cref="Analysis.ForceVector" />
		private Vector<double>? _currentForces;

		/// <summary>
		///     The residual force <see cref="Vector" /> of current iteration
		/// </summary>
		/// <inheritdoc cref="Analysis.ForceVector" />
		private Vector<double>? _currentResidual;

		/// <summary>
		///     The secant stiffness <see cref="Matrix" /> of current iteration
		/// </summary>
		private Matrix<double>? _currentStiffness;

		/// <summary>
		///     Field to store current iteration.
		/// </summary>
		private int _iteration;

		/// <summary>
		///     The displacement <see cref="Vector" /> of last iteration
		/// </summary>
		/// <inheritdoc cref="Analysis.DisplacementVector" />
		private Vector<double>? _lastDisplacements;

		/// <summary>
		///     The residual <see cref="Vector" /> of last iteration
		/// </summary>
		/// <inheritdoc cref="Analysis.ForceVector" />
		private Vector<double>? _lastResidual;

		/// <summary>
		///     The secant stiffness <see cref="Matrix" /> of last iteration
		/// </summary>
		private Matrix<double>? _lastStiffness;

		/// <summary>
		///     Field to store current load step.
		/// </summary>
		private int _loadStep;

		/// <summary>
		///     Monitored displacements list.
		/// </summary>
		/// <inheritdoc cref="Analysis.DisplacementVector" />
		private List<MonitoredDisplacement>? _monitoredDisplacements;

		/// <summary>
		///     Field to store the DoF index for <see cref="_monitoredDisplacements" />.
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

				for (var i = 0; i < _currentResidual!.Count; i++)
				{
					num += _currentResidual[i] * _currentResidual[i];
					den += _currentForces![i] * _currentForces[i];
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
		///     Get current load factor.
		/// </summary>
		private double LoadFactor => (double) _loadStep / NumLoadSteps;

		#endregion

		#region Constructors

		/// <summary>
		///     Secant analysis constructor.
		/// </summary>
		/// <param name="numLoadSteps">The number of load steps to perform (default: 50).</param>
		/// <param name="tolerance">The convergence tolerance (default: 1E-6).</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step (default: 10000).</param>
		/// <param name="minIterations">Minimum number of iterations for each load step (default: 2).</param>
		/// <inheritdoc />
		public SecantAnalysis(FEMInput femInput, int numLoadSteps = 50, double tolerance = 1E-6, int maxIterations = 10000, int minIterations = 2)
			: base(femInput)
		{
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
		/// <param name="monitoredIndex">The index of a degree of freedom to monitor, if wanted.</param>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis.ForceVector" /> (default: 1).</param>
		public void Execute(int? monitoredIndex = null, double loadFactor = 1)
		{
			// Get force vector
			ForceVector = FemInput.ForceVector * loadFactor;

			// Get the initial stiffness and force vector simplified
			UpdateStiffness();

			// Initiate lists
			Initiate(monitoredIndex);

			// Analysis by load steps
			StepAnalysis();

			// Set displacements
			DisplacementVector = _currentDisplacements;
			GlobalStiffness    = _currentStiffness;

			// Set Reactions
			FemInput.Grips.SetReactions(GetReactions());
		}


		/// <summary>
		///     Generate an <see cref="FEMOutput" /> from analysis results.
		/// </summary>
		/// <returns>
		///     null if no monitored index was provided.
		/// </returns>
		public FEMOutput? GenerateOutput() => _monitoredDisplacements is null
			? null
			: new FEMOutput(_monitoredDisplacements);

		/// <summary>
		///     Update displacements.
		/// </summary>
		private void DisplacementUpdate()
		{
			// Set last displacements
			_lastDisplacements = _currentDisplacements!.Clone();

			// Increment displacements
			_currentDisplacements += _currentStiffness!.Solve(-_currentResidual);

			// Update displacements in grips
			FemInput.Grips.SetDisplacements(_currentDisplacements);
		}

		/// <summary>
		///     Initiate fields.
		/// </summary>
		/// <inheritdoc cref="Execute" />
		private void Initiate(int? monitoredIndex)
		{
			_monitoredIndex = monitoredIndex;

			if (_monitoredIndex.HasValue)
				_monitoredDisplacements = new List<MonitoredDisplacement>();

			// Calculate initial stiffness
			_currentStiffness = GlobalStiffness;

			// Calculate initial displacements
			var lf0 = (double) 1 / NumLoadSteps;
			_currentDisplacements = _currentStiffness!.Solve(lf0 * ForceVector);

			// Initiate solution values
			_lastStiffness     = _currentStiffness.Clone();
			_lastDisplacements = Vector<double>.Build.Dense(FemInput.NumberOfDoFs);
			_lastResidual      = Vector<double>.Build.Dense(FemInput.NumberOfDoFs);
			_currentResidual   = Vector<double>.Build.Dense(FemInput.NumberOfDoFs);
		}

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		private void Iterate()
		{
			// Initiate first iteration
			_iteration = 1;

			while (_iteration <= MaxIterations)
			{
				// Calculate element forces
				FemInput.Elements.CalculateForces();

				// Update residual
				ResidualUpdate();

				// Check convergence or if analysis must stop
				if (ConvergenceReached || StopCheck())
					return;

				// Update stiffness and displacements
				SecantStiffnessUpdate();
				DisplacementUpdate();

				// Increase iteration count
				_iteration++;
			}
		}

		/// <summary>
		///     Calculate residual force <see cref="Vector" />.
		/// </summary>
		private Vector<double> ResidualForces() => InternalForces(FemInput) - _currentForces;

		/// <summary>
		///     Update residual force <see cref="Vector" />.
		/// </summary>
		private void ResidualUpdate()
		{
			// Set new values
			_lastResidual    = _currentResidual!.Clone();
			_currentResidual = ResidualForces();
		}

		/// <summary>
		///     Save load step results after achieving convergence.
		/// </summary>
		private void SaveLoadStepResults()
		{
			if (!_monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = Length.FromMillimeters(_currentDisplacements![_monitoredIndex.Value]);

			// Add to list
			_monitoredDisplacements!.Add(new MonitoredDisplacement(disp, LoadFactor));
		}

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix" /> of current iteration.
		/// </summary>
		/// <param name="simplify">Simplify stiffness?</param>
		private void SecantStiffnessUpdate(bool simplify = true)
		{
			// Clone current stiffness
			var kCur = _currentStiffness!.Clone();

			// Calculate the variation of displacements and residual as vectors
			Vector<double>
				dDisp = _currentDisplacements - _lastDisplacements,
				dRes  = _currentResidual - _lastResidual;

			// Increment current stiffness
			var dK = ((dRes - _lastStiffness * dDisp) / dDisp.Norm(2)).ToColumnMatrix() * dDisp.ToRowMatrix();

			// Set new values
			_currentStiffness = _lastStiffness + dK;
			_lastStiffness    = kCur;

			// Simplify
			if (simplify)
				Simplify(simplify, false);
		}

		/// <summary>
		///     Execute step by step analysis.
		/// </summary>
		private void StepAnalysis()
		{
			// Initiate first load step
			_loadStep = 1;

			while (_loadStep <= NumLoadSteps)
			{
				// Get the force vector
				_currentForces = LoadFactor * ForceVector;

				// Iterate
				Iterate();

				// Verify if convergence was not reached
				if (Stop)
					return;

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
			Stop = _iteration == MaxIterations || _currentResidual!.ContainsNaN() ||
			       _currentDisplacements!.ContainsNaN() || _currentStiffness!.ContainsNaN();

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
		private bool VerifyConvergence(double convergence) => convergence <= Tolerance && _iteration >= MinIterations;

		#endregion

	}
}