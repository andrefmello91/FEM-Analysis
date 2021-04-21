using System;
using andrefmello91.OnPlaneComponents;
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
		/// <summary>
		///		Linear-elastic analysis.
		/// </summary>
		Linear,
		
		/// <summary>
		///		Nonlinear analysis.
		/// </summary>
		Nonlinear
	}

	/// <summary>
	///		Nonlinear solution procedures.
	/// </summary>
	public enum NonLinearSolver
	{
		/// <summary>
		///		Newton-Raphson nonlinear solver.
		/// </summary>
		NewtonRaphson,
		
		/// <summary>
		///		Modified Newton-Raphson nonlinear solver.
		/// </summary>
		ModifiedNewtonRaphson,
		
		/// <summary>
		///		Secant Method nonlinear solver.
		/// </summary>
		Secant
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
		/// <remarks>
		///     Components in <see cref="LengthUnit.Millimeter" />.
		/// </remarks>
		Vector<double> Displacements { get; set; }

		/// <summary>
		///     Get the force <see cref="Vector" /> in this element, in global coordinate system.
		/// </summary>
		/// <remarks>
		///     Components in <see cref="ForceUnit.Newton" />.
		/// </remarks>
		Vector<double> Forces { get; set; }

		/// <summary>
		///     Get the grips of this element.
		/// </summary>
		IGrip[] Grips { get; }

		/// <summary>
		///     Get stiffness <see cref="Matrix" /> in the global coordinate system.
		/// </summary>
		/// <remarks>
		///     Components in units compatible to <see cref="ForceUnit.Newton" /> and <see cref="LengthUnit.Millimeter" />.
		/// </remarks>
		Matrix<double> Stiffness { get; set; }

		#endregion

		#region Methods

		/// <summary>
		///     Calculate forces in this element after updating displacements at each grip.
		/// </summary>
		void CalculateForces();

		#endregion

	}

	/// <summary>
	///		Interface for nonlinear finite elements.
	/// </summary>
	public interface INonlinearElement : IFiniteElement
	{
		/// <summary>
		///		The results of the last iteration.
		/// </summary>
		IterationResult LastIterationResult { get; set; }
		
		/// <summary>
		///		The results of the current iteration.
		/// </summary>
		IterationResult CurrentIterationResult { get; set; }
	}
}