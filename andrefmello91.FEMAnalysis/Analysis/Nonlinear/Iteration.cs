using System.Collections.Generic;
using andrefmello91.Extensions;
using andrefmello91.FEMAnalysis;
using MathNet.Numerics.LinearAlgebra;
using static andrefmello91.FEMAnalysis.NonlinearAnalysis;

namespace andrefmello91.FEMAnalysis
{

	/// <summary>
	///     Class for nonlinear iteration results.
	/// </summary>
	public class Iteration : IIteration, ICloneable<Iteration>
	{
		#region Properties

		/// <inheritdoc/>
		public double DisplacementConvergence { get; protected set; }

		/// <inheritdoc/>
		public virtual Vector<double> DisplacementIncrement { get; private set; }

		/// <inheritdoc/>
		public Vector<double> Displacements { get; protected set; }

		/// <inheritdoc/>
		public double ForceConvergence { get; protected set; }

		/// <inheritdoc/>
		public Vector<double> InternalForces { get; set; }

		/// <inheritdoc/>
		public int Number { get; set; }

		/// <inheritdoc/>
		public Vector<double> ResidualForces { get; private set; }

		/// <inheritdoc/>
		public Matrix<double> Stiffness { get; set; }

		#endregion

		#region Constructors

		/// <inheritdoc cref="From(int,bool)"/>
		protected Iteration(int numberOfDoFs)
			: this(Vector<double>.Build.Dense(numberOfDoFs), Vector<double>.Build.Dense(numberOfDoFs), Matrix<double>.Build.Dense(numberOfDoFs, numberOfDoFs))
		{
		}

		/// <inheritdoc cref="From(Vector{double}, Vector{double}, Matrix{double}, bool)"/>
		protected Iteration(Vector<double> displacements, Vector<double> residualForces, Matrix<double> stiffness)
		{
			Displacements         = displacements;
			ResidualForces        = residualForces;
			Stiffness             = stiffness;
			InternalForces        = Vector<double>.Build.Dense(displacements.Count);
			DisplacementIncrement = Vector<double>.Build.Dense(displacements.Count);
		}

		#endregion

		#region Methods

		/// <summary>
		///     Create an iteration object.
		/// </summary>
		/// <param name="numberOfDoFs">The number of degrees of freedom.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public static IIteration From(int numberOfDoFs, bool simulate = false) => simulate switch
		{
			false => new Iteration(numberOfDoFs),
			_     => new SimulationIteration(numberOfDoFs)
		};
		
		/// <summary>
		///     Create an iteration object.
		/// </summary>
		/// <param name="displacements">The displacement vector of this iteration.</param>
		/// <param name="residualForces">The residual force vector of this iteration.</param>
		/// <param name="stiffness">The stiffness matrix of this iteration.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public static IIteration From(Vector<double> displacements, Vector<double> residualForces, Matrix<double> stiffness, bool simulate = false) => simulate switch
		{
			false => new Iteration(displacements, residualForces, stiffness),
			_     => new SimulationIteration(displacements, residualForces, stiffness)
		};
		
		/// <summary>
		///     Create an iteration from a load step result.
		/// </summary>
		/// <param name="loadStep">A calculated load step.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public static IIteration FromStepResult(LoadStep loadStep, bool simulate = false) => simulate switch
		{
			false => new Iteration(loadStep.FinalDisplacements, Vector<double>.Build.Dense(loadStep.FinalDisplacements.Count), loadStep.Stiffness),
			_     => new SimulationIteration(loadStep.FinalDisplacements, Vector<double>.Build.Dense(loadStep.FinalDisplacements.Count), loadStep.Stiffness)
		};

		/// <inheritdoc/>
		public void CalculateConvergence(IEnumerable<double> appliedForces, IEnumerable<double> initialIncrement)
		{
			ForceConvergence        = NonlinearAnalysis.CalculateConvergence(ResidualForces, appliedForces);
			DisplacementConvergence = NonlinearAnalysis.CalculateConvergence(DisplacementIncrement, initialIncrement);
		}

		/// <inheritdoc/>
		public virtual bool CheckConvergence(AnalysisParameters parameters) =>
			Number >= parameters.MinIterations &&
			ForceConvergence <= parameters.ForceTolerance &&
			DisplacementConvergence <= parameters.DisplacementTolerance;

		/// <inheritdoc/>
		public bool CheckStopCondition(AnalysisParameters parameters) =>
			Number >= parameters.MaxIterations || ResidualForces.ContainsNaNOrInfinity() ||
			Displacements.ContainsNaNOrInfinity() || Stiffness.ContainsNaN();

		/// <inheritdoc/>
		public void IncrementDisplacements(Vector<double> displacementIncrement)
		{
			DisplacementIncrement =  displacementIncrement;
			Displacements         += displacementIncrement;
		}

		/// <inheritdoc/>
		public void UpdateForces(Vector<double> appliedForces, Vector<double> internalForces)
		{
			InternalForces = internalForces;
			ResidualForces = internalForces - appliedForces;
		}

		#region Interface Implementations

		/// <inheritdoc />
		public Iteration Clone() => new(Displacements.Clone(), ResidualForces.Clone(), Stiffness.Clone()) { Number = Number };

		#endregion

		#region Object override

		/// <inheritdoc />
		IIteration ICloneable<IIteration>.Clone() => Clone();

		/// <inheritdoc />
		public override string ToString() => $"Iteration {Number}";

		#endregion

		#endregion

		#region Operators

		/// <summary>
		///     Get the number of a iteration.
		/// </summary>
		/// <returns>
		///     <see cref="Iteration.Number" />
		/// </returns>
		public static explicit operator int(Iteration iteration) => iteration.Number;

		/// <summary>
		///     Check the iteration number.
		/// </summary>
		/// <returns>
		///     True if the iteration number is equal to the right number.
		/// </returns>
		public static bool operator ==(Iteration left, int right) => left.Number == right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is bigger than the right number.
		/// </returns>
		public static bool operator >(Iteration left, int right) => left.Number > right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is smaller than the right number.
		/// </returns>
		public static bool operator <(Iteration left, int right) => left.Number < right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is bigger or equal to the right number.
		/// </returns>
		public static bool operator >=(Iteration left, int right) => left.Number >= right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is smaller or equal to the right number.
		/// </returns>
		public static bool operator <=(Iteration left, int right) => left.Number <= right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is not equal to the right number.
		/// </returns>
		public static bool operator !=(Iteration left, int right) => left.Number != right;

		#endregion

	}
}