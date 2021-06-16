using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet.Units;
#nullable disable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Finite element input class.
	/// </summary>
	/// <typeparam name="TFiniteElement">Any type that implements <see cref="IFiniteElement" />.</typeparam>
	public interface IFEMInput<out TFiniteElement> : IEnumerable<TFiniteElement>
		where TFiniteElement : IFiniteElement
	{

		#region Properties

		/// <summary>
		///     Get the index of constrained degrees of freedom.
		/// </summary>
		public List<int> ConstraintIndex { get; }

		/// <summary>
		///     Get the external force <see cref="Vector" />.
		/// </summary>
		/// <remarks>
		///     Components in <see cref="ForceUnit.Newton" />.
		/// </remarks>
		public Vector<double> ForceVector { get; }

		/// <summary>
		///     Get the grips of the finite element model.
		/// </summary>
		public List<IGrip> Grips { get; }

		/// <summary>
		///     Get the number of degrees of freedom (DoFs).
		/// </summary>
		public int NumberOfDoFs { get; }

		#endregion

	}

	/// <summary>
	///     Default finite element input class.
	/// </summary>
	/// <inheritdoc cref="IFEMInput{TFiniteElement}" />
	public class FEMInput<TFiniteElement> : List<TFiniteElement>, IFEMInput<TFiniteElement>
		where TFiniteElement : IFiniteElement
	{

		#region Properties

		#region Interface Implementations

		/// <inheritdoc />
		public List<int> ConstraintIndex { get; }

		/// <inheritdoc />
		public Vector<double> ForceVector { get; }

		/// <inheritdoc />
		public List<IGrip> Grips { get; }

		/// <inheritdoc />
		public int NumberOfDoFs { get; }

		#endregion

		#endregion

		#region Constructors

		/// <inheritdoc cref="FEMInput{IFiniteElement}(IEnumerable{IFiniteElement}, IEnumerable{IGrip})" />
		/// <remarks>
		///     Grips are taken from <paramref name="elements" />.
		/// </remarks>
		public FEMInput(IEnumerable<TFiniteElement> elements)
			: this(elements, elements.SelectMany(e => e.Grips).Distinct().OrderBy(g => g.Number).ToList())
		{
		}

		/// <summary>
		///     Input Data constructor.
		/// </summary>
		/// <param name="elements">The collection containing all distinct <see cref="IFiniteElement" />'s in the model.</param>
		/// <param name="grips">The collection containing all distinct <see cref="IGrip" />'s in the model.</param>
		public FEMInput(IEnumerable<TFiniteElement> elements, IEnumerable<IGrip> grips)
			: base(elements)
		{
			Grips           = grips.ToList();
			NumberOfDoFs    = 2 * Grips.Count;
			ForceVector     = Grips.AssembleForceVector();
			ConstraintIndex = Grips.GetConstraintIndex().ToList();
		}

		#endregion

		#region Methods

		#region Object override

		/// <inheritdoc />
		public override string ToString() =>
			$"Number of elements: {Count}\n" +
			$"Number of grips: {Grips.Count}\n" +
			$"Force vector: \n{ForceVector}\n" +
			$"Constraint Index: {ConstraintIndex.Select(i => i.ToString()).Aggregate((i, f) => $"{i} - {f}")}";

		#endregion

		#endregion

	}
}