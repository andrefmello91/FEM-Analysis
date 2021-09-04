using System;
using andrefmello91.OnPlaneComponents;
using MathNet.Numerics.LinearAlgebra.Double;
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

		/// <inheritdoc cref="IFEMInput.Forces" />
		/// <remarks>
		///     Simplified by constrained DoFs.
		/// </remarks>
		public ForceVector Forces { get; protected set; }

		/// <summary>
		///     Get/set global stiffness <see cref="Matrix" />.
		/// </summary>
		public StiffnessMatrix GlobalStiffness { get; protected set; }

		#endregion

		/// <summary>
		///     Event to execute when analysis is complete.
		/// </summary>
		public abstract event EventHandler? AnalysisComplete;

		/// <summary>
		///     Event to execute when analysis is aborted.
		/// </summary>
		public abstract event EventHandler? AnalysisAborted;

		#region Constructors

		/// <summary>
		///     Base analysis constructor.
		/// </summary>
		/// <param name="femInput">The <see cref="IFEMInput" /> for finite element analysis.</param>
		protected Analysis(IFEMInput femInput)
		{
			FemInput        = femInput;
			Displacements   = DisplacementVector.Zero(FemInput.NumberOfDoFs);
			Forces          = FemInput.Forces.Simplified(FemInput.ConstraintIndex);
			GlobalStiffness = femInput.AssembleStiffness();
		}

		#endregion

		#region Methods

		/// <summary>
		///     Calculate the <see cref="Vector" /> of support reactions.
		/// </summary>
		public ForceVector GetReactions() => (ForceVector) (GlobalStiffness * Displacements - Forces);

		/// <summary>
		///     Invoke the event.
		/// </summary>
		/// <param name="handler">The handler.</param>
		protected void Invoke(EventHandler? handler) => handler?.Invoke(this, EventArgs.Empty);

		/// <inheritdoc cref="Invoke" />
		protected void Invoke<TEventArgs>(EventHandler<TEventArgs>? handler, TEventArgs? eventArgs) where TEventArgs : EventArgs =>
			handler.Invoke(this, eventArgs);

		/// <inheritdoc />
		public override string ToString() =>
			$"{FemInput}\n" +
			"Global Stiffness:\n" +
			$"{GlobalStiffness}\n" +
			"Displacement Vector:\n" +
			$"{Displacements}";

		#endregion

	}
}