using System;

namespace andrefmello91.FEMAnalysis;

/// <summary>
///     Event args for load steps.
/// </summary>
public class StepEventArgs : EventArgs
{

	/// <summary>
	///     The current load step.
	/// </summary>
	public LoadStep Step { get; }

	/// <summary>
	///     Create an event arg for load step.
	/// </summary>
	/// <param name="step">Current load step</param>
	public StepEventArgs(LoadStep step) => Step = step;
}