using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;

namespace andrefmello91.FEMAnalysis
{

	/// <summary>
	///     Class for nonlinear iteration results.
	/// </summary>
	public class Iteration : IIteration, ICloneable<Iteration>
	{

		#region Properties

		/// <inheritdoc />
		public double DisplacementConvergence { get; protected set; }

		/// <inheritdoc />
		public virtual DisplacementVector DisplacementIncrement { get; private set; }

		/// <inheritdoc />
		public DisplacementVector Displacements { get; protected set; }

		/// <inheritdoc />
		public double ForceConvergence { get; protected set; }

		/// <inheritdoc />
		public ForceVector InternalForces { get; set; }

		/// <inheritdoc />
		public int Number { get; set; }

		/// <inheritdoc />
		public ForceVector ResidualForces { get; private set; }

		/// <inheritdoc />
		public StiffnessMatrix Stiffness { get; set; }

		#endregion

		#region Constructors

		/// <inheritdoc cref="From(int,bool)" />
		protected Iteration(int numberOfDoFs)
			: this(DisplacementVector.Zero(numberOfDoFs), ForceVector.Zero(numberOfDoFs), StiffnessMatrix.Zero(numberOfDoFs))
		{
		}

		/// <inheritdoc cref="From(DisplacementVector, ForceVector, StiffnessMatrix, bool)" />
		protected Iteration(DisplacementVector displacements, ForceVector residualForces, StiffnessMatrix stiffness)
		{
			Displacements         = displacements;
			ResidualForces        = residualForces;
			Stiffness             = stiffness;
			InternalForces        = ForceVector.Zero(displacements.Count);
			DisplacementIncrement = DisplacementVector.Zero(displacements.Count);
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
		/// <param name="loadFactor">The load factor if <paramref name="simulate" /> is true.</param>
		public static IIteration From(DisplacementVector displacements, ForceVector residualForces, StiffnessMatrix stiffness, bool simulate = false, double loadFactor = 0) => simulate switch
		{
			false => new Iteration(displacements, residualForces, stiffness),
			_     => new SimulationIteration(displacements, residualForces, stiffness, loadFactor)
		};

		/// <summary>
		///     Create an iteration from a load step result.
		/// </summary>
		/// <param name="loadStep">A calculated load step.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public static IIteration FromStepResult(LoadStep loadStep, bool simulate = false) => simulate switch
		{
			false => new Iteration(loadStep.FinalDisplacements, ForceVector.Zero(loadStep.FinalDisplacements.Count), loadStep.Stiffness),
			_     => new SimulationIteration(loadStep.FinalDisplacements, ForceVector.Zero(loadStep.FinalDisplacements.Count), loadStep.Stiffness, loadStep.LoadFactor)
		};

		/// <inheritdoc />
		public override string ToString() => $"Iteration {Number}";

		/// <inheritdoc />
		public Iteration Clone() => new((DisplacementVector) Displacements.Clone(), (ForceVector) ResidualForces.Clone(), (StiffnessMatrix) Stiffness.Clone()) { Number = Number };

		/// <inheritdoc />
		public void CalculateConvergence(ForceVector appliedForces, DisplacementVector initialIncrement)
		{
			ForceConvergence        = NonlinearAnalysis.CalculateConvergence(ResidualForces, appliedForces);
			DisplacementConvergence = NonlinearAnalysis.CalculateConvergence(DisplacementIncrement, initialIncrement);
		}

		/// <inheritdoc />
		public virtual bool CheckConvergence(AnalysisParameters parameters) =>
			Number >= parameters.MinIterations &&
			ForceConvergence <= parameters.ForceTolerance && DisplacementConvergence <= parameters.DisplacementTolerance;

		/// <inheritdoc />
		public bool CheckStopCondition(AnalysisParameters parameters) =>
			Number >= parameters.MaxIterations || ResidualForces.ContainsNaNOrInfinity() ||
			Displacements.ContainsNaNOrInfinity() || Stiffness.ContainsNaN();

		/// <inheritdoc />
		public void IncrementDisplacements(DisplacementVector displacementIncrement)
		{
			DisplacementIncrement = displacementIncrement;
			Displacements         = (DisplacementVector) (Displacements + displacementIncrement);
		}

		/// <inheritdoc />
		public void UpdateForces(ForceVector appliedForces, ForceVector internalForces)
		{
			InternalForces = internalForces;
			ResidualForces = (ForceVector) (internalForces - appliedForces);
		}

		/// <inheritdoc />
		IIteration ICloneable<IIteration>.Clone() => Clone();

		#endregion

		#region Operators

		/// <summary>
		///     Check the iteration number.
		/// </summary>
		/// <returns>
		///     True if the iteration number is equal to the right number.
		/// </returns>
		public static bool operator ==(Iteration left, int right) => left.Number == right;

		/// <summary>
		///     Get the number of a iteration.
		/// </summary>
		/// <returns>
		///     <see cref="Iteration.Number" />
		/// </returns>
		public static explicit operator int(Iteration iteration) => iteration.Number;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is bigger than the right number.
		/// </returns>
		public static bool operator >(Iteration left, int right) => left.Number > right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is bigger or equal to the right number.
		/// </returns>
		public static bool operator >=(Iteration left, int right) => left.Number >= right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is not equal to the right number.
		/// </returns>
		public static bool operator !=(Iteration left, int right) => left.Number != right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is smaller than the right number.
		/// </returns>
		public static bool operator <(Iteration left, int right) => left.Number < right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the iteration number is smaller or equal to the right number.
		/// </returns>
		public static bool operator <=(Iteration left, int right) => left.Number <= right;

		#endregion

	}
}