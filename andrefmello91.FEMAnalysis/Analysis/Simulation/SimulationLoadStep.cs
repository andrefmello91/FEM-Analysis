using System;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.FEMAnalysis.Simulation;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;
using static andrefmello91.FEMAnalysis.Analysis<andrefmello91.FEMAnalysis.IFiniteElement>;
using static andrefmello91.FEMAnalysis.NonlinearAnalysis;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Class for load step results.
	/// </summary>
	public class SimulationLoadStep : LoadStep<SimulationIteration>, ICloneable<SimulationLoadStep>
	{

		#region Constructors

		/// <summary>
		///     Create a load step object.
		/// </summary>
		/// <param name="number">The number of this step.</param>
		/// <param name="numberOfDoFs">The number of degrees of freedom.</param>
		private SimulationLoadStep(int numberOfDoFs, AnalysisParameters parameters, int number = 0)
			: this(number, Vector<double>.Build.Dense(numberOfDoFs), Vector<double>.Build.Dense(numberOfDoFs), Matrix<double>.Build.Dense(numberOfDoFs, numberOfDoFs), parameters)
		{
		}

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="loadFactor">The load factor of this step.</param>
		/// <param name="number">The number of this step.</param>
		/// <param name="forces">The force vector of this step.</param>
		private SimulationLoadStep(Vector<double> forces, double loadFactor, AnalysisParameters parameters, int number = 0)
			: this(number, forces, Vector<double>.Build.Dense(forces.Count), Matrix<double>.Build.Dense(forces.Count, forces.Count), parameters) =>
			LoadFactor = loadFactor;

		/// <inheritdoc cref="LoadStep" />
		/// <param name="initialDisplacements">The initial displacement vector of this step.</param>
		/// <param name="stiffness">The stiffness matrix of this step.</param>
		/// <param name="parameters">The analysis parameters.</param>
		private SimulationLoadStep(int number, Vector<double> forces, Vector<double> initialDisplacements, Matrix<double> stiffness, AnalysisParameters parameters)
			: base(number, forces, initialDisplacements, stiffness, parameters)
		{
			Iterations.Add(new SimulationIteration(initialDisplacements, Vector<double>.Build.Dense(initialDisplacements.Count), stiffness));
		}

		#endregion

		#region Methods

		/// <summary>
		///     Create a load step for a simulation procedure.
		/// </summary>
		/// <inheritdoc cref="LoadStep.From"/>
		public static SimulationLoadStep From(IFEMInput<IFiniteElement> femInput, double loadFactor, AnalysisParameters parameters, int stepNumber) =>
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
		public static SimulationLoadStep InitialStep(IFEMInput<IFiniteElement> femInput, double loadFactor, AnalysisParameters parameters)
		{
			var step      = From(femInput, loadFactor, parameters, 1);
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

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		/// <param name="femInput">The finite element input.</param>
		public virtual void Iterate(IFEMInput<IFiniteElement> femInput)
		{
			// Initiate first iteration
			foreach (var iteration in this)
				iteration.Number = 0;

			// Iterate
			do
			{
				// Add iteration
				NewIteration();

				// Update stiffness and displacements
				UpdateDisplacements(this, femInput);
				UpdateStiffness(this, femInput);

				// Calculate element forces
				femInput.CalculateForces();

				// Update internal forces
				var extForces = SimplifiedForces(Forces, femInput.ConstraintIndex);
				var intForces = femInput.AssembleInternalForces();
				OngoingIteration.UpdateForces(extForces, intForces);

				// Calculate convergence
				OngoingIteration.CalculateConvergence(extForces, FirstIteration.DisplacementIncrement);
				
			} while (!IterativeStop());
		}

		#region Interface Implementations

		/// <inheritdoc />
		public SimulationLoadStep Clone() => new(Number, Forces.Clone(), FinalDisplacements.Clone(), Stiffness.Clone(), Parameters)
		{
			LoadFactor = LoadFactor
		};

		#endregion

		#region Object override

		/// <inheritdoc />
		public override string ToString() => $"Load step {Number}";

		#endregion

		#endregion
		
	}
}