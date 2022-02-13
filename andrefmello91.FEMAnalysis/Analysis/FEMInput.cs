using System.Collections.Generic;
using System.Linq;
using andrefmello91.OnPlaneComponents;
using MathNet.Numerics.LinearAlgebra.Double;
#nullable disable

namespace andrefmello91.FEMAnalysis;

/// <summary>
///     Finite element input class.
/// </summary>
public interface IFEMInput : IEnumerable<IFiniteElement>
{

	/// <summary>
	///     Get the index of constrained degrees of freedom.
	/// </summary>
	List<int> ConstraintIndex { get; }

	/// <summary>
	///     Get the external force <see cref="Vector" />.
	/// </summary>
	ForceVector Forces { get; }

	/// <summary>
	///     Get the grips of the finite element model.
	/// </summary>
	List<IGrip> Grips { get; }

	/// <summary>
	///     The monitored elements.
	/// </summary>
	IEnumerable<IMonitoredElement> MonitoredElements { get; }

	/// <summary>
	///     Get the number of degrees of freedom (DoFs).
	/// </summary>
	int NumberOfDoFs { get; }
}

/// <summary>
///     Default finite element input class.
/// </summary>
/// <inheritdoc cref="IFEMInput{TFiniteElement}" />
public class FEMInput : List<IFiniteElement>, IFEMInput
{

	/// <inheritdoc />
	public List<int> ConstraintIndex { get; }

	/// <inheritdoc />
	public ForceVector Forces { get; }

	/// <inheritdoc />
	public List<IGrip> Grips { get; }

	/// <summary>
	///     The monitored elements.
	/// </summary>
	public IEnumerable<IMonitoredElement> MonitoredElements => this
		.Where(e => e is IMonitoredElement { Monitored: true })
		.Cast<IMonitoredElement>()
		.Concat(Grips
			.Where(e => e is IMonitoredElement { Monitored: true })
			.Cast<IMonitoredElement>());

	/// <inheritdoc />
	public int NumberOfDoFs { get; }

	/// <inheritdoc cref="FEMInput(IEnumerable{IFiniteElement}, IEnumerable{IGrip})" />
	/// <remarks>
	///     Grips are taken from <paramref name="elements" />.
	/// </remarks>
	public FEMInput(IEnumerable<IFiniteElement> elements)
		: this(elements, elements.SelectMany(e => e.Grips).Distinct())
	{
	}

	/// <summary>
	///     Input Data constructor.
	/// </summary>
	/// <param name="elements">The collection containing all distinct <see cref="IFiniteElement" />'s in the model.</param>
	/// <param name="grips">The collection containing all distinct <see cref="IGrip" />'s in the model.</param>
	public FEMInput(IEnumerable<IFiniteElement> elements, IEnumerable<IGrip> grips)
		: base(elements)
	{
		Grips           = grips.OrderBy(g => g.Number).ToList();
		NumberOfDoFs    = 2 * Grips.Count;
		ConstraintIndex = Grips.GetConstraintIndex().ToList();
		Forces          = this.AssembleExternalForces(false);
	}

	/// <inheritdoc />
	public override string ToString() =>
		$"Number of elements: {Count}\n" +
		$"Number of grips: {Grips.Count}\n" +
		$"Force vector: \n{Forces}\n" +
		$"Constraint Index: {ConstraintIndex.Select(i => i.ToString()).Aggregate((i, f) => $"{i} - {f}")}";
}