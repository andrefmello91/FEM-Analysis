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
		///		The monitored displacement of this load step.
		/// </summary>
		public MonitoredDisplacement? MonitoredDisplacement { get; set; }
		
		/// <summary>
		///		The displacement vector of this load step.
		/// </summary>
		public Vector<double> Displacements { get; set; }
		
		/// <summary>
		///		The force vector of this load step.
		/// </summary>
		public Vector<double> Forces { get; set; }
		
		/// <summary>
		///		The stiffness matrix of this load step.
		/// </summary>
		public Matrix<double> Stiffness { get; set; }

		///  <summary>
		/// 		Create a load step object.
		///  </summary>
		///  <param name="number">The number of this load step.</param>
		///  <param name="forces">The force vector of this load step.</param>
		///  <param name="displacements">The displacement vector of this load step.</param>
		///  <param name="stiffness">The stiffness matrix of this load step.</param>
		public LoadStepResult(int number, Vector<double> forces, Vector<double> displacements, Matrix<double> stiffness)
		{
			Number        = number;
			Forces        = forces;
			Displacements = displacements;
			Stiffness     = stiffness;
		}

		/// <inheritdoc />
		public LoadStepResult Clone() => new(Number, Forces.Clone(), Displacements.Clone(), Stiffness.Clone());

		/// <summary>
		///		Get the number of a load step.
		/// </summary>
		/// <returns>
		///		<see cref="LoadStepResult.Number"/>
		/// </returns>
		public static explicit operator int(LoadStepResult loadStepResult) => loadStepResult.Number;
	}
}