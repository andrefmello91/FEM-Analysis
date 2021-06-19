using System.Collections.Generic;
using andrefmello91.Extensions;
using andrefmello91.FEMAnalysis.Simulation;
using MathNet.Numerics.LinearAlgebra;
using static andrefmello91.FEMAnalysis.NonlinearAnalysis;

namespace andrefmello91.FEMAnalysis
{
	public interface IIteration : ICloneable<IIteration>
	{
		/// <summary>
		///     The displacement convergence of this iteration.
		/// </summary>
		double DisplacementConvergence { get; }

		/// <summary>
		///     The displacement increment vector from external forces of this iteration.
		/// </summary>
		Vector<double> DisplacementIncrement { get; }

		/// <summary>
		///     The displacement vector of this iteration.
		/// </summary>
		Vector<double> Displacements { get; }

		/// <summary>
		///     The force convergence of this iteration.
		/// </summary>
		double ForceConvergence { get; }

		/// <summary>
		///     The internal force vector of this iteration.
		/// </summary>
		Vector<double> InternalForces { get; }

		/// <summary>
		///     The number of this iteration.
		/// </summary>
		int Number { get; set; }

		/// <summary>
		///     The residual force vector of this iteration.
		/// </summary>
		Vector<double> ResidualForces { get; }

		/// <summary>
		///     The stiffness matrix of this iteration.
		/// </summary>
		Matrix<double> Stiffness { get; set; }

		/// <summary>
		///     Calculate the convergence of this iteration.
		/// </summary>
		/// <param name="appliedForces">The applied forces of the current step.</param>
		/// <param name="initialIncrement">The displacement increment of the first iteration.</param>
		void CalculateConvergence(IEnumerable<double> appliedForces, IEnumerable<double> initialIncrement);

		/// <summary>
		///     Check convergence for this iteration.
		/// </summary>
		/// <param name="parameters">The analysis parameters.</param>
		/// <returns>
		///     True if this iteration number is equal or bigger than minimum iterations and force or displacement convergences are
		///     smaller than their respective tolerances.
		/// </returns>
		bool CheckConvergence(AnalysisParameters parameters);

		/// <summary>
		///     Check the stop condition for this iteration.
		/// </summary>
		/// <inheritdoc cref="CheckConvergence" />
		/// <returns>
		///     True if this iteration number is equal or bigger than maximum number of iterations or any of tha analysis vectors
		///     and matrix contains <see cref="double.NaN" />.
		/// </returns>
		bool CheckStopCondition(AnalysisParameters parameters);

		/// <summary>
		///     Increment displacements of this iteration.
		/// </summary>
		/// <param name="displacementIncrement">The vector of displacement increments.</param>
		void IncrementDisplacements(Vector<double> displacementIncrement);

		/// <summary>
		///     Update forces in this iteration.
		/// </summary>
		/// <param name="appliedForces">The vector of applied forces of the current step.</param>
		/// <param name="internalForces">The vector of internal forces.</param>
		void UpdateForces(Vector<double> appliedForces, Vector<double> internalForces);
	}
	
	/// <summary>
	///     Class for nonlinear iteration results.
	/// </summary>
	public class Iteration : IIteration, ICloneable<Iteration>
	{

		#region Properties

		/// <inheritdoc/>
		public double DisplacementConvergence { get; private set; }

		/// <inheritdoc/>
		public Vector<double> DisplacementIncrement { get; private set; }

		/// <inheritdoc/>
		public Vector<double> Displacements { get; protected set; }

		/// <inheritdoc/>
		public double ForceConvergence { get; private set; }

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

		/// <summary>
		///     Create an iteration object with elements composed by zero..
		/// </summary>
		/// <param name="numberOfDoFs">The number of degrees of freedom.</param>
		public Iteration(int numberOfDoFs)
			: this(Vector<double>.Build.Dense(numberOfDoFs), Vector<double>.Build.Dense(numberOfDoFs), Matrix<double>.Build.Dense(numberOfDoFs, numberOfDoFs))
		{
		}

		/// <summary>
		///     Create an iteration object.
		/// </summary>
		/// <param name="displacements">The displacement vector of this iteration.</param>
		/// <param name="residualForces">The residual force vector of this iteration.</param>
		/// <param name="stiffness">The stiffness matrix of this iteration.</param>
		public Iteration(Vector<double> displacements, Vector<double> residualForces, Matrix<double> stiffness)
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
		///     Create an iteration result from a load step result.
		/// </summary>
		/// <param name="loadStep">A calculated load step.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public static IIteration FromStepResult<TIteration>(LoadStep<TIteration> loadStep, bool simulate = false)
			where TIteration : class, IIteration => simulate switch
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
		public bool CheckConvergence(AnalysisParameters parameters) =>
			Number >= parameters.MinIterations &&
			(ForceConvergence <= parameters.ForceTolerance || DisplacementConvergence <= parameters.DisplacementTolerance);

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