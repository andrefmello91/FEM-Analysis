using System;
using andrefmello91.OnPlaneComponents;
using andrefmello91.OnPlaneComponents.Displacement;
using andrefmello91.OnPlaneComponents.Force;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     The analysis types.
	/// </summary>
	public enum AnalysisType
	{
		Linear,
		NonLinear
	}

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
		///     Get the displacement <see cref="Vector" />, in global coordinate system.
		/// </summary>
		/// <inheritdoc cref="LocalDisplacements" />
		Vector<double> Displacements { get; }

		/// <summary>
		///     Get the force <see cref="Vector" /> in this element, in global coordinate system.
		/// </summary>
		/// <inheritdoc cref="LocalForces" />
		Vector<double> Forces { get; }

		/// <summary>
		///     Get the grips of this element.
		/// </summary>
		IGrip[] Grips { get; }

		/// <summary>
		///     Get the displacement <see cref="Vector" />, in local coordinate system.
		/// </summary>
		/// <remarks>
		///     Components in <see cref="LengthUnit.Millimeter" />.
		/// </remarks>
		Vector<double> LocalDisplacements { get; }

		/// <summary>
		///     Get the force <see cref="Vector" />, in local coordinate system.
		/// </summary>
		/// <remarks>
		///     Components in <see cref="ForceUnit.Newton" />.
		/// </remarks>
		Vector<double> LocalForces { get; }

		/// <summary>
		///     Get local stiffness <see cref="Matrix" />.
		/// </summary>
		Matrix<double> LocalStiffness { get; }

		/// <summary>
		///     Get global stiffness <see cref="Matrix" />.
		/// </summary>
		Matrix<double> Stiffness { get; }

		/// <summary>
		///     Get the transformation <see cref="Matrix" />.
		/// </summary>
		Matrix<double> TransformationMatrix { get; }

		#endregion
		#region Methods

		/// <summary>
		///     Calculate forces in this element after updating displacements at each grip.
		/// </summary>
		void CalculateForces();

		#endregion
	}
}