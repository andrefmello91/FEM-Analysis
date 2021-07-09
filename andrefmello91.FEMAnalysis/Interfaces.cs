using System;
using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;

namespace andrefmello91.FEMAnalysis
{

	/// <summary>
	///     Interface for numbered elements.
	/// </summary>
	public interface INumberedElement
	{

		#region Properties

		/// <summary>
		///     Get the index of degrees of freedom of this element.
		/// </summary>
		int[] DoFIndex { get; }

		/// <summary>
		///     Get/set the number of this element.
		/// </summary>
		/// <remarks>
		///     Enumeration starts at 1.
		/// </remarks>
		int Number { get; set; }

		#endregion

	}

	/// <summary>
	///     Interface for grips.
	/// </summary>
	public interface IGrip : INumberedElement, IEquatable<IGrip>, IComparable<IGrip>
	{

		#region Properties

		/// <summary>
		///     Get the <see cref="Constraint" /> in this grip.
		/// </summary>
		Constraint Constraint { get; }

		/// <summary>
		///     Get the <see cref="PlaneDisplacement" /> in this grip.
		/// </summary>
		PlaneDisplacement Displacement { get; set; }

		/// <summary>
		///     Get the <see cref="PlaneForce" /> in this grip.
		/// </summary>
		PlaneForce Force { get; }

		/// <summary>
		///     Get the <see cref="PlaneForce" /> reaction in this grip, in case it's constrained.
		/// </summary>
		PlaneForce Reaction { get; set; }

		#endregion

	}

	/// <summary>
	///     Interface for finite elements.
	/// </summary>
	public interface IFiniteElement : INumberedElement, IEquatable<IFiniteElement>, IComparable<IFiniteElement>
	{

		#region Properties

		/// <summary>
		///     Get the displacement vector, in global coordinate system.
		/// </summary>
		DisplacementVector Displacements { get; set; }

		/// <summary>
		///     Get the force vector in this element, in global coordinate system.
		/// </summary>
		ForceVector Forces { get; set; }

		/// <summary>
		///     Get the grips of this element.
		/// </summary>
		IGrip[] Grips { get; }

		/// <summary>
		///     Get stiffness matrix in the global coordinate system.
		/// </summary>
		StiffnessMatrix Stiffness { get; set; }

		#endregion

		#region Methods

		/// <summary>
		///     Calculate forces in this element after updating displacements at each grip.
		/// </summary>
		void CalculateForces();

		/// <summary>
		///     Update displacements in this element.
		/// </summary>
		void UpdateDisplacements();

		/// <summary>
		///     Update stiffness of this element.
		/// </summary>
		void UpdateStiffness();

		#endregion

	}

	/// <summary>
	///     Interface for iterations.
	/// </summary>
	public interface IIteration : ICloneable<IIteration>
	{

		#region Properties

		/// <summary>
		///     The displacement convergence of this iteration.
		/// </summary>
		double DisplacementConvergence { get; }

		/// <summary>
		///     The displacement increment vector from external forces of this iteration.
		/// </summary>
		DisplacementVector DisplacementIncrement { get; }

		/// <summary>
		///     The displacement vector of this iteration.
		/// </summary>
		DisplacementVector Displacements { get; }

		/// <summary>
		///     The force convergence of this iteration.
		/// </summary>
		double ForceConvergence { get; }

		/// <summary>
		///     The internal force vector of this iteration.
		/// </summary>
		ForceVector InternalForces { get; }

		/// <summary>
		///     The number of this iteration.
		/// </summary>
		int Number { get; set; }

		/// <summary>
		///     The residual force vector of this iteration.
		/// </summary>
		ForceVector ResidualForces { get; }

		/// <summary>
		///     The stiffness matrix of this iteration.
		/// </summary>
		StiffnessMatrix Stiffness { get; set; }

		#endregion

		#region Methods

		/// <summary>
		///     Calculate the convergence of this iteration.
		/// </summary>
		/// <param name="appliedForces">The applied forces of the current step.</param>
		/// <param name="initialIncrement">The displacement increment of the first iteration.</param>
		void CalculateConvergence(ForceVector appliedForces, DisplacementVector initialIncrement);

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
		void IncrementDisplacements(DisplacementVector displacementIncrement);

		/// <summary>
		///     Update forces in this iteration.
		/// </summary>
		/// <param name="appliedForces">The vector of applied forces of the current step.</param>
		/// <param name="internalForces">The vector of internal forces.</param>
		void UpdateForces(ForceVector appliedForces, ForceVector internalForces);

		#endregion

	}
}