using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Nonlinear analysis class
	/// </summary>
	public class NonlinearAnalysis : Analysis
	{
		/// <summary>
		///		The results of the last iteration.
		/// </summary>
		private IterationResult _lastResult;

		/// <summary>
		///		The results of the current iteration.
		/// </summary>
		private IterationResult _currentResult;

		/// <summary>
		///		The list of load step results.
		/// </summary>
		private List<LoadStepResult> _loadSteps = new();
		
		/// <summary>
		///     The force vector of current load step.
		/// </summary>
		/// <inheritdoc cref="Analysis.ForceVector" />
		private Vector<double> _currentForces;

		/// <summary>
		///     Monitored displacements list.
		/// </summary>
		/// <inheritdoc cref="Analysis.DisplacementVector" />
		private List<MonitoredDisplacement>? _monitoredDisplacements;

		/// <summary>
		///     Field to store the DoF index for <see cref="_monitoredDisplacements" />.
		/// </summary>
		private int? _monitoredIndex;

		/// <summary>
		///     Current iteration number.
		/// </summary>
		private int _iteration;
		
		/// <summary>
		///     Current load step number.
		/// </summary>
		private int _loadStep;

		/// <summary>
		///		The nonlinear equation solver.
		/// </summary>
		public NonLinearSolver Solver { get; set; }

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

				var res = _currentResult.ResidualForces;
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
		///     Get current load factor.
		/// </summary>
		private double LoadFactor => (double) _loadStep / NumLoadSteps;
		
		/// <summary>
		///		Nonlinear analysis constructor.
		/// </summary>
		/// <param name="nonlinearInput">The <see cref="NonlinearInput"/>.</param>
		/// <param name="solver">The <see cref="NonLinearSolver"/> to use.</param>
		/// <param name="numLoadSteps">The number of load steps to perform (default: 50).</param>
		/// <param name="tolerance">The convergence tolerance (default: 1E-6).</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step (default: 10000).</param>
		/// <param name="minIterations">Minimum number of iterations for each load step (default: 2).</param>
		public NonlinearAnalysis(
			NonlinearInput nonlinearInput,
			NonLinearSolver solver = NonLinearSolver.NewtonRaphson,
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
		
		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="convergence">
		///     Calculated convergence.
		///     <para>See: <see cref="Convergence" />.</para>
		/// </param>
		private bool VerifyConvergence(double convergence) => convergence <= Tolerance && _iteration >= MinIterations;

	}
}