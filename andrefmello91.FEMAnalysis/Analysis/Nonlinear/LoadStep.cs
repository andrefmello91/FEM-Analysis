﻿using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;

using static andrefmello91.FEMAnalysis.Analysis<andrefmello91.FEMAnalysis.IFiniteElement>;

namespace andrefmello91.FEMAnalysis
{


	/// <summary>
	///     Class for load step results.
	/// </summary>
	public class LoadStep : LoadStep<Iteration>, ICloneable<LoadStep>
	{
		#region Constructors

		/// <summary>
		///     Create a load step object.
		/// </summary>
		/// <param name="number">The number of this step.</param>
		/// <param name="numberOfDoFs">The number of degrees of freedom.</param>
		private LoadStep(int numberOfDoFs, AnalysisParameters parameters, int number = 0)
			: this(number, Vector<double>.Build.Dense(numberOfDoFs), Vector<double>.Build.Dense(numberOfDoFs), Matrix<double>.Build.Dense(numberOfDoFs, numberOfDoFs), parameters)
		{
		}

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="loadFactor">The load factor of this step.</param>
		/// <param name="number">The number of this step.</param>
		/// <param name="forces">The force vector of this step.</param>
		private LoadStep(Vector<double> forces, double loadFactor, AnalysisParameters parameters, int number = 0)
			: this(number, forces, Vector<double>.Build.Dense(forces.Count), Matrix<double>.Build.Dense(forces.Count, forces.Count), parameters) =>
			LoadFactor = loadFactor;

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="initialDisplacements">The initial displacement vector of this step.</param>
		/// <param name="stiffness">The stiffness matrix of this step.</param>
		/// <param name="parameters">The analysis parameters.</param>
		private LoadStep(int number, Vector<double> forces, Vector<double> initialDisplacements, Matrix<double> stiffness, AnalysisParameters parameters)
			: base(number, forces, initialDisplacements, stiffness, parameters)
		{
			Iterations.Add(new Iteration(initialDisplacements, Vector<double>.Build.Dense(initialDisplacements.Count), stiffness));
		}

		#endregion

		#region Methods

		/// <summary>
		///     Create a load step for a nonlinear analysis procedure.
		/// </summary>
		/// <param name="femInput">The finite element input.</param>
		/// <param name="loadFactor">The load factor of the first step. </param>
		/// <param name="parameters">The analysis parameters.</param>
		/// <param name="stepNumber">The number of the load step.</param>
		public static LoadStep From(IFEMInput<IFiniteElement> femInput, double loadFactor, AnalysisParameters parameters, int stepNumber) =>
			new(femInput.ForceVector * loadFactor, loadFactor, parameters, stepNumber);

		/// <summary>
		///     Do the initial load step of a nonlinear analysis procedure.
		/// </summary>
		/// <param name="femInput">The finite element input.</param>
		/// <param name="loadFactor">The load factor of the first step. </param>
		/// <param name="parameters">The analysis parameters.</param>
		/// <returns>
		///     The initial <see cref="LoadStep" />.
		/// </returns>
		public static LoadStep InitialStep(IFEMInput<IFiniteElement> femInput, double loadFactor, AnalysisParameters parameters)
		{
			var step = From(femInput, loadFactor, parameters, 1);
			
			var iteration = step.OngoingIteration;

			// Get the initial stiffness and force vector simplified
			iteration.Stiffness = femInput.AssembleStiffness();
			var stiffness = SimplifiedStiffness(iteration.Stiffness, femInput.ConstraintIndex);

			// Calculate initial displacements
			var fi = SimplifiedForces(step.Forces, femInput.ConstraintIndex);
			iteration.IncrementDisplacements(stiffness.Solve(fi));

			// Update displacements in grips and elements
			femInput.Grips.SetDisplacements(iteration.Displacements);
			femInput.UpdateDisplacements();

			// Calculate element forces
			femInput.CalculateForces();

			// Update internal forces
			iteration.InternalForces = femInput.AssembleInternalForces();

			return step;
		}


		#region Interface Implementations

		/// <inheritdoc />
		public LoadStep Clone() => new(Number, Forces.Clone(), FinalDisplacements.Clone(), Stiffness.Clone(), Parameters)
		{
			LoadFactor = LoadFactor
		};

		#endregion
		
		#endregion

	}
}