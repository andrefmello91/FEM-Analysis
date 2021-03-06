using System;
using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;

namespace andrefmello91.FEMAnalysis;

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
	///     The name of this element.
	/// </summary>
	string Name { get; }

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
	DisplacementVector Displacements { get; }

	/// <summary>
	///     Get the force vector in this element, in global coordinate system.
	/// </summary>
	ForceVector Forces { get; }

	/// <summary>
	///     Get the grips of this element.
	/// </summary>
	IGrip[] Grips { get; }

	/// <summary>
	///     Get stiffness matrix in the global coordinate system.
	/// </summary>
	StiffnessMatrix Stiffness { get; }

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
///     Interface for monitored elements.
/// </summary>
public interface IMonitoredElement
{

	#region Properties

	/// <summary>
	///     The element monitor.
	/// </summary>
	ElementMonitor? Monitor { get; }

	/// <summary>
	///     Get/set monitoring of this element.
	/// </summary>
	public bool Monitored { get; set; }

	#endregion

	#region Methods

	/// <summary>
	///     Add the monitored values to monitor for the current load factor.
	/// </summary>
	/// <param name="loadFactor">The current load factor.</param>
	void AddValue(double loadFactor);

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
	/// <inheritdoc cref="LoadStep.FullForceVector" />
	ForceVector InternalForces { get; }

	/// <summary>
	///     The number of this iteration.
	/// </summary>
	int Number { get; set; }

	/// <summary>
	///     The residual force vector of this iteration.
	/// </summary>
	/// <inheritdoc cref="LoadStep.FullForceVector" />
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
	/// <param name="appliedForces">The applied forces of the current step, simplified in constrained DoFs.</param>
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
	/// <remarks>
	///     Forces must be simplified in constrained DoFs.
	/// </remarks>
	/// <param name="appliedForces">The vector of applied forces of the current step, simplified in constrained DoFs.</param>
	/// <param name="internalForces">The vector of internal forces, simplified in constrained DoFs.</param>
	void UpdateForces(ForceVector appliedForces, ForceVector internalForces);

	#endregion

}

/// <summary>
///     Generic interface for monitored values.
/// </summary>
/// <typeparam name="TStruct">The monitored type.</typeparam>
public interface IMonitoredValue<TStruct> : IEquatable<IMonitoredValue<TStruct>>, IComparable<IMonitoredValue<TStruct>>
	where TStruct : struct
{

	#region Properties

	/// <summary>
	///     The load factor associated to value.
	/// </summary>
	double LoadFactor { get; }

	/// <summary>
	///     The monitored value.
	/// </summary>
	TStruct Value { get; }

	#endregion

}

public interface IElementMonitor
{
}