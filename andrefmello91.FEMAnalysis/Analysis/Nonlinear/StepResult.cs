using System;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.FEMAnalysis.Simulation;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Class for step results.
	/// </summary>
	public class StepResult : List<IterationResult>, ICloneable<StepResult>
	{

		#region Properties

		/// <summary>
		///     The convergence of this step.
		/// </summary>
		public double Convergence { get; set; }
		
		/// <summary>
		///     The load factor of this step.
		/// </summary>
		public double LoadFactor { get; private set; }

		/// <summary>
		///     The displacement vector of this step.
		/// </summary>
		public Vector<double> Displacements { get; private set; }

		/// <summary>
		///     The force vector of this step.
		/// </summary>
		public Vector<double> Forces { get; private set; }

		/// <summary>
		///     The status of this step. True if it was calculated.
		/// </summary>
		public bool IsCalculated { get; set; }

		/// <summary>
		///     The monitored displacement of this step.
		/// </summary>
		public MonitoredDisplacement? MonitoredDisplacement { get; set; }

		/// <summary>
		///     The number of this step.
		/// </summary>
		public int Number { get; set; }

		/// <summary>
		///     The stiffness matrix of this step.
		/// </summary>
		public Matrix<double> Stiffness { get; private set; }
		
		#endregion

		#region Constructors

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="number">The number of this step.</param>
		/// <param name="numberOfDoFs">The number of degrees of freedom.</param>
		public StepResult(int numberOfDoFs, int number = 0)
			: this(number, Vector<double>.Build.Dense(numberOfDoFs), Vector<double>.Build.Dense(numberOfDoFs), Matrix<double>.Build.Dense(numberOfDoFs, numberOfDoFs))
		{
		}

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="loadFactor">The load factor of this step.</param>
		/// <param name="number">The number of this step.</param>
		/// <param name="forces">The force vector of this step.</param>
		public StepResult(Vector<double> forces, double loadFactor, int number = 0)
			: this(number, forces, Vector<double>.Build.Dense(forces.Count), Matrix<double>.Build.Dense(forces.Count, forces.Count))
		{
			LoadFactor = loadFactor;
		}

		/// <inheritdoc cref="StepResult" />
		/// <param name="displacements">The displacement vector of this step.</param>
		/// <param name="stiffness">The stiffness matrix of this step.</param>
		public StepResult(int number, Vector<double> forces, Vector<double> displacements, Matrix<double> stiffness)
		{
			Number        = number;
			Forces        = forces;
			Displacements = displacements;
			Stiffness     = stiffness;
		}

		#endregion

		#region Methods

		#region Interface Implementations

		/// <summary>
		///		Add a new iteration in this load step.
		/// </summary>
		///  <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public void NewIteration(bool simulate = false) => Add(this.Any() ? this.Last().Clone() : IterationResult.FromStepResult(this, simulate));
		
		/// <summary>
		///		Increment forces in this step.
		/// </summary>
		/// <param name="loadFactorIncrement">The increment of the load factor.</param>
		public void IncrementLoad(double loadFactorIncrement)
		{
			if (this.Any() && this.Last() is SimulationIterationResult itResult)
				itResult.LoadFactorIncrement = loadFactorIncrement;
			
			// Get the actual force multiplier
			var lf = 1D + loadFactorIncrement / LoadFactor;
			
			// Update values
			LoadFactor += loadFactorIncrement;
			Forces     *= lf;
		}
		
		/// <summary>
		///     Set step results after achieving convergence.
		/// </summary>
		public void SetResults(int? monitoredIndex = null)
		{
			IsCalculated  = true;
			Convergence   = this.Last().ForceConvergence;
			Displacements = this.Last().Displacements;
			Stiffness     = this.Last().Stiffness;

			if (!monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = Length.FromMillimeters(Displacements[monitoredIndex.Value]);

			// Set to step
			MonitoredDisplacement = new MonitoredDisplacement(disp, LoadFactor);
		}

		/// <summary>
		///		Get the accumulated displacement increment at this load step.
		/// </summary>
		/// <param name="finalIndex">The final index to consider increments.</param>
		public Vector<double> AccumulatedDisplacementIncrement(Index finalIndex) =>
			this[finalIndex].Displacements - Find(i => i == 1)!.Displacements;

		/// <inheritdoc />
		public StepResult Clone() => new(Number, Forces.Clone(), Displacements.Clone(), Stiffness.Clone())
		{
			LoadFactor = LoadFactor
		};

		/// <inheritdoc />
		public override string ToString() => $"Step {Number}";

		#endregion

		#endregion

		#region Operators

		/// <summary>
		///     Get the number of a step.
		/// </summary>
		/// <returns>
		///     <see cref="StepResult.Number" />
		/// </returns>
		public static explicit operator int(StepResult stepResult) => stepResult.Number;

		/// <summary>
		///		Check the step number.
		/// </summary>
		/// <returns>
		///		True if the step number is equal to the right number.
		/// </returns>
		public static bool operator ==(StepResult left, int right) => left.Number == right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is not equal to the right number.
		/// </returns>
		public static bool operator !=(StepResult left, int right) => left.Number != right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is bigger than the right number.
		/// </returns>
		public static bool operator >(StepResult left, int right) => left.Number > right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is smaller than the right number.
		/// </returns>
		public static bool operator <(StepResult left, int right) => left.Number < right;
		
		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is bigger or equal to the right number.
		/// </returns>
		public static bool operator >=(StepResult left, int right) => left.Number >= right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is smaller or equal to the right number.
		/// </returns>
		public static bool operator <=(StepResult left, int right) => left.Number <= right;

		#endregion

	}
}