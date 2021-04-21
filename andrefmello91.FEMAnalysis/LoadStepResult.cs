using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Class for load step results.
	/// </summary>
	public class LoadStepResult : ICloneable<LoadStepResult>
	{
		/// <summary>
		///		The number of this load step.
		/// </summary>
		public int Number { get; set; }
		
		/// <summary>
		///		The displacement vector of this load step.
		/// </summary>
		public Vector<double>? Displacements { get; set; }
		
		/// <summary>
		///		The force vector of this load step.
		/// </summary>
		public Vector<double> Forces { get; set; }

		/// <summary>
		///		Create a load step object.
		/// </summary>
		/// <param name="forces">The force vector of this load step.</param>
		public LoadStepResult(Vector<double> forces)
		{
			Forces = forces;
		}

		/// <inheritdoc />
		public LoadStepResult Clone() => new(Forces.Clone())
		{
			Displacements = Displacements?.Clone()
		};

		/// <summary>
		///		Get the number of a load step.
		/// </summary>
		/// <returns>
		///		<see cref="LoadStepResult.Number"/>
		/// </returns>
		public static explicit operator int(LoadStepResult loadStepResult) => loadStepResult.Number;
	}
}