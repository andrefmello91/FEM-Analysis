using System.Collections.Generic;
using Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

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
		private Vector<double> _currentDisplacements;

		/// <summary>
		///     Field to store each iteration force <see cref="Vector" />
		/// </summary>
		/// <inheritdoc cref="Analysis.ForceVector" />
		private Vector<double> _currentForces;

		/// <summary>
		///     The residual force <see cref="Vector" /> of current iteration
		/// </summary>
		/// <inheritdoc cref="Analysis.ForceVector" />
		private Vector<double> _currentResidual;

		/// <summary>
		///     The secant stiffness <see cref="Matrix" /> of current iteration
		/// </summary>
		private Matrix<double> _currentStiffness;

		/// <summary>
		///     Field to store current iteration.
		/// </summary>
		private int _iteration;

		/// <summary>
		///     The displacement <see cref="Vector" /> of last iteration
		/// </summary>
		/// <inheritdoc cref="Analysis.DisplacementVector" />
		private Vector<double> _lastDisplacements;

		/// <summary>
		///     The residual <see cref="Vector" /> of last iteration
		/// </summary>
		/// <inheritdoc cref="Analysis.ForceVector" />
		private Vector<double> _lastResidual;

		/// <summary>
		///     The secant stiffness <see cref="Matrix" /> of last iteration
		/// </summary>
		private Matrix<double> _lastStiffness;

		/// <summary>
		///     Field to store current load step.
		/// </summary>
		private int _loadStep;

		/// <summary>
		///     Field to store the maximum number of iterations.
		/// </summary>
		private int _maxIterations;

		/// <summary>
		///     Field to store the minimum number of iterations.
		/// </summary>
		private int _minIterations;

		/// <summary>
		///     Field to store the DoF index for <see cref="MonitoredDisplacements" />.
		/// </summary>
		private int? _monitoredIndex;

		/// <summary>
		///     Field to store the number of load steps to perform.
		/// </summary>
		private int _numLoadSteps;

		/// <summary>
		///     Field to store convergence tolerance.
		/// </summary>
		private double _tolerance;

		#endregion

		#region Properties

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

				for (var i = 0; i < _currentResidual.Count; i++)
				{
					num += _currentResidual[i] * _currentResidual[i];
					den += _currentForces[i] * _currentForces[i];
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
		private double LoadFactor => (double) _loadStep / _numLoadSteps;

		/// <summary>
		///     Get/set values of displacements monitored.
		/// </summary>
		/// <inheritdoc cref="Analysis.DisplacementVector" />
		public List<double> MonitoredDisplacements { get; private set; }

		/// <summary>
		///     Get/set values of load factor associated to <see cref="MonitoredDisplacements" />.
		/// </summary>
		public List<double> MonitoredLoadFactor { get; private set; }

		/// <summary>
		///     Get/set when to stop analysis.
		/// </summary>
		public bool Stop { get; private set; }

		/// <summary>
		///     Get/set the stop message.
		/// </summary>
		public string StopMessage { get; private set; }

		#endregion

		#region Constructors

		/// <summary>
		///     Secant analysis object
		/// </summary>
		/// <inheritdoc />
		public SecantAnalysis(InputData inputData)
			: base(inputData)
		{
		}

		#endregion

		#region  Methods

		/// <summary>
		///     Do the analysis.
		/// </summary>
		/// <param name="monitoredIndex">The DoF index to monitor, if wanted.</param>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis.ForceVector" /> (default: 1).</param>
		/// <param name="numLoadSteps">The number of load steps to perform (default: 50).</param>
		/// <param name="tolerance">The convergence tolerance (default: 1E-6).</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step (default: 10000).</param>
		/// <param name="minIterations">Minimum number of iterations for each load step (default: 2).</param>
		public void Do(int? monitoredIndex = null, double loadFactor = 1, int numLoadSteps = 50, double tolerance = 1E-6, int maxIterations = 10000, int minIterations = 2)
		{
			// Get force vector
			ForceVector = InputData.ForceVector * loadFactor;

			// Get the initial stiffness and force vector simplified
			UpdateStiffness();

			// Initiate lists
			Initiate(monitoredIndex, numLoadSteps, tolerance, maxIterations, minIterations);

			// Analysis by load steps
			StepAnalysis();

			// Set displacements
			DisplacementVector = _currentDisplacements;
			GlobalStiffness = _currentStiffness;
			NodalDisplacements(_currentDisplacements);
		}

		/// <summary>
		///     Initiate fields.
		/// </summary>
		/// <param name="monitoredIndex">The DoF index to monitor, if wanted.</param>
		/// <param name="numLoadSteps">The number of load steps to perform.</param>
		/// <param name="tolerance">The convergence tolerance.</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step.</param>
		/// <param name="minIterations">Minimum number of iterations for each load step.</param>
		private void Initiate(int? monitoredIndex, int numLoadSteps, double tolerance, int maxIterations, int minIterations)
		{
			_monitoredIndex = monitoredIndex;

			MonitoredDisplacements = _monitoredIndex.HasValue ? new List<double>() : null;
			MonitoredLoadFactor    = _monitoredIndex.HasValue ? new List<double>() : null;

			_numLoadSteps  = numLoadSteps;
			_tolerance     = tolerance;
			_maxIterations = maxIterations;
			_minIterations = minIterations;

			// Calculate initial stiffness
			_currentStiffness = GlobalStiffness;

			// Calculate initial displacements
			var lf0 = (double) 1 / _numLoadSteps;
			_currentDisplacements = _currentStiffness.Solve(lf0 * ForceVector);

			// Initiate solution values
			_lastStiffness     = _currentStiffness.Clone();
			_lastDisplacements = Vector<double>.Build.Dense(NumberOfDoFs);
			_lastResidual      = Vector<double>.Build.Dense(NumberOfDoFs);
			_currentResidual   = Vector<double>.Build.Dense(NumberOfDoFs);
		}

		/// <summary>
		///     Do analysis by load steps.
		/// </summary>
		private void StepAnalysis()
		{
			for (_loadStep = 1; _loadStep <= _numLoadSteps; _loadStep++)
			{
				// Get the force vector
				_currentForces = LoadFactor * ForceVector;

				// Iterate
				Iterate();

				// Verify if convergence was not reached
				if (Stop)
					break;

				// Set load step results
				SaveLoadStepResults();
			}
		}

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		private void Iterate()
		{
			for (_iteration = 1; _iteration <= _maxIterations; _iteration++)
			{
				// Calculate element displacements and forces
				ElementAnalysis(_currentDisplacements);

				// Update residual
				ResidualUpdate();

				// Check convergence or if analysis must stop
				if (ConvergenceReached || StopCheck())
					break;

				// Update stiffness and displacements
				SecantStiffnessUpdate();
				DisplacementUpdate();
			}
		}

		/// <summary>
		///     Check if analysis must stop.
		/// </summary>
		private bool StopCheck()
		{
			// Check if one stop condition is reached
			Stop = _iteration == _maxIterations || _currentResidual.ContainsNaN() || _currentDisplacements.ContainsNaN() || _currentStiffness.ContainsNaN();

			// Check if maximum number of iterations is reached
			if (Stop)
				StopMessage = $"Convergence not reached at load step {_loadStep}";

			return Stop;
		}

		/// <summary>
		///     Update residual force <see cref="Vector" />.
		/// </summary>
		private void ResidualUpdate()
		{
			// Set new values
			_lastResidual    = _currentResidual.Clone();
			_currentResidual = ResidualForces();
		}

		/// <summary>
		///     Calculate residual force <see cref="Vector" />.
		/// </summary>
		private Vector<double> ResidualForces() => InternalForces() - _currentForces;

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix" /> of current iteration.
		/// </summary>
		/// <param name="simplify">Simplify stiffness?</param>
		private void SecantStiffnessUpdate(bool simplify = true)
		{
			// Clone current stiffness
			var kCur = _currentStiffness.Clone();

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
				SimplifyStiffness();
		}

		/// <summary>
		///     Simplify the secant stiffness <see cref="Matrix" /> of current iteration.
		/// </summary>
		private void SimplifyStiffness()
		{
			_currentStiffness.ClearRows(ConstraintIndex);
			_currentStiffness.ClearColumns(ConstraintIndex);

			foreach (var i in ConstraintIndex)
				_currentStiffness[i, i] = 1;
		}

		/// <summary>
		///     Update displacements.
		/// </summary>
		private void DisplacementUpdate()
		{
			// Set last displacements
			_lastDisplacements = _currentDisplacements.Clone();

			// Increment displacements
			_currentDisplacements += _currentStiffness.Solve(-_currentResidual);
		}


		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="convergence">Calculated convergence.
		///     <para>See: <see cref="Convergence" />.</para>
		/// </param>
		private bool VerifyConvergence(double convergence) => convergence <= _tolerance && _iteration >= _minIterations;

		/// <summary>
		///     Save load step results after achieving convergence.
		/// </summary>
		private void SaveLoadStepResults()
		{
			if (!_monitoredIndex.HasValue)
				return;

			MonitoredDisplacements.Add(_currentDisplacements[_monitoredIndex.Value]);
			MonitoredLoadFactor.Add(LoadFactor);
		}

		#endregion
	}
}