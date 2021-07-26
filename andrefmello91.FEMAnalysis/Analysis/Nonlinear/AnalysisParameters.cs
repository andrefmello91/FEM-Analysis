namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Nonlinear analysis parameters struct.
	/// </summary>
	/// <param name="Solver">The nonlinear equation solver.</param>
	/// <param name="NumberOfSteps">The number of steps to execute.</param>
	/// <param name="MaxIterations">The maximum number of iterations.</param>
	/// <param name="MinIterations">The minimum number of iterations.</param>
	/// <param name="ForceTolerance">The convergence tolerance for residual forces.</param>
	/// <param name="DisplacementTolerance">The convergence tolerance for displacement increments.</param>
	public record AnalysisParameters(NonLinearSolver Solver, int NumberOfSteps, int MaxIterations, int MinIterations, double ForceTolerance, double DisplacementTolerance)
	{

		#region Properties

		/// <summary>
		///     Get the default analysis parameters.
		/// </summary>
		/// <remarks>
		///     Solver: <see cref="NonLinearSolver.NewtonRaphson" />.
		///     <para>NumberOfSteps: 50.</para>
		///     <para>MaxIterations: 10000.</para>
		///     <para>MinIterations: 2.</para>
		///     <para>ForceTolerance: 1E-3</para>
		///     <para>DisplacementTolerance: 1E-8</para>
		/// </remarks>
		public static AnalysisParameters Default { get; } = new(NonLinearSolver.NewtonRaphson, 50, 1000, 2, 1E-3, 1E-8);

		/// <summary>
		///     Get/set the convergence tolerance for displacement increments.
		/// </summary>
		/// <remarks>
		///     Default: 1E-6
		/// </remarks>
		public double DisplacementTolerance { get; set; } = DisplacementTolerance;

		/// <summary>
		///     Get/set the convergence tolerance for residual forces.
		/// </summary>
		/// <remarks>
		///     Default: 1E-3
		/// </remarks>
		public double ForceTolerance { get; set; } = ForceTolerance;

		/// <summary>
		///     Get/set the maximum number of iterations.
		/// </summary>
		/// <remarks>
		///     Default: 1000
		/// </remarks>
		public int MaxIterations { get; set; } = MaxIterations;

		/// <summary>
		///     Get/set the minimum number of iterations.
		/// </summary>
		/// <remarks>
		///     Default: 2
		/// </remarks>
		public int MinIterations { get; set; } = MinIterations;

		/// <summary>
		///     Get/set the number of steps to execute.
		/// </summary>
		/// <remarks>
		///     Default: 50
		/// </remarks>
		public int NumberOfSteps { get; set; } = NumberOfSteps;

		/// <summary>
		///     The nonlinear equation solver.
		/// </summary>
		public NonLinearSolver Solver { get; set; } = Solver;

		#endregion

	}
}