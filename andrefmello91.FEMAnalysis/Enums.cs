namespace andrefmello91.FEMAnalysis;

/// <summary>
///     The analysis types.
/// </summary>
public enum AnalysisType
{
	/// <summary>
	///     Linear-elastic analysis.
	/// </summary>
	Linear,

	/// <summary>
	///     Nonlinear analysis.
	/// </summary>
	Nonlinear
}

/// <summary>
///     Nonlinear solution procedures.
/// </summary>
public enum NonLinearSolver
{
	/// <summary>
	///     Newton-Raphson nonlinear solver.
	/// </summary>
	NewtonRaphson,

	/// <summary>
	///     Modified Newton-Raphson nonlinear solver.
	/// </summary>
	ModifiedNewtonRaphson,

	/// <summary>
	///     Secant Method nonlinear solver.
	/// </summary>
	Secant
}

/// <summary>
///     Types of nonlinear analysis control.
/// </summary>
public enum AnalysisControl
{
	/// <summary>
	///     Control the analysis by increasing the applied forces.
	/// </summary>
	/// <remarks>
	///     This cannot predict the full response of the model.
	/// </remarks>
	Force,

	/// <summary>
	///     Control the analysis by increasing the displacements.
	/// </summary>
	/// <remarks>
	///     This can predict the full response of the model.
	/// </remarks>
	Displacement
}