using System.Collections.Generic;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet;
using UnitsNet.Units;

#nullable enable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Analysis base class.
	/// </summary>
	public abstract class Analysis
	{

		#region Properties

		/// <summary>
		///     Get/set the displacement <see cref="Vector" />.
		/// </summary>
		/// <remarks>
		///     Components in <see cref="LengthUnit.Millimeter" />.
		/// </remarks>
		public DisplacementVector Displacements { get; protected set; }

		/// <summary>
		///     Get the input for finite element analysis.
		/// </summary>
		public IFEMInput FemInput { get; }

		/// <inheritdoc cref="IFEMInput.ForceVector" />
		public ForceVector Forces { get; protected set; }

		/// <summary>
		///     Get/set global stiffness <see cref="Matrix" />.
		/// </summary>
		public StiffnessMatrix GlobalStiffness { get; protected set; }

		#endregion

		#region Constructors

		/// <summary>
		///     Base analysis constructor.
		/// </summary>
		/// <param name="femInput">The <see cref="IFEMInput" /> for finite element analysis.</param>
		protected Analysis(IFEMInput femInput)
		{
			FemInput                      = femInput;
			Displacements                 = DisplacementVector.Zero(FemInput.NumberOfDoFs);
			Forces                        = FemInput.ForceVector;
			GlobalStiffness               = StiffnessMatrix.Assemble(FemInput);
			Displacements.ConstraintIndex = FemInput.ConstraintIndex;
		}

		#endregion

		#region Methods

		/// <summary>
		///     Calculate the <see cref="Vector" /> of support reactions.
		/// </summary>
		public ForceVector GetReactions() => GlobalStiffness * Displacements - Forces;

		#region Object override

		/// <inheritdoc />
		public override string ToString() =>
			$"{FemInput}\n" +
			"Global Stiffness:\n" +
			$"{GlobalStiffness}\n" +
			"Displacement Vector:\n" +
			$"{Displacements}";

		#endregion

		#endregion

	}
}