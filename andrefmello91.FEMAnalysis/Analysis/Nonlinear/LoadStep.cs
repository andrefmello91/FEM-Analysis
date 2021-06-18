using System;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.FEMAnalysis.Simulation;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

using static andrefmello91.FEMAnalysis.Analysis<andrefmello91.FEMAnalysis.IFiniteElement>;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Class for load step results.
	/// </summary>
	public class LoadStep : List<IterationResult>, ICloneable<LoadStep>
	{

		#region Properties

		/// <summary>
		///		The analysis parameters.
		/// </summary>
		public AnalysisParameters Parameters { get; }

		/// <summary>
		///     The convergence of this step.
		/// </summary>
		public double Convergence { get; set; }
		
		/// <summary>
		///     The load factor of this step.
		/// </summary>
		public double LoadFactor { get; private set; }

		/// <summary>
		///     The displacement vector at the beginning of this step.
		/// </summary>
		public Vector<double> InitialDisplacements { get; }

		/// <summary>
		///     The displacement vector at the end of this step.
		/// </summary>
		public Vector<double> FinalDisplacements => this.Last().Displacements;

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
		///     The current stiffness matrix of this step (stiffness of the current iteration.
		/// </summary>
		public Matrix<double> Stiffness => this.Last().Stiffness;
		
		/// <summary>
		///		Get the first iteration of the current step.
		/// </summary>
		public IterationResult FirstIteration => Find(i => (int) i == 1)!;

		/// <summary>
		///     The results of the ongoing iteration.
		/// </summary>
		public IterationResult OngoingIteration => this[^1];

		/// <summary>
		///     The results of the current solution (last solved iteration [i - 1]).
		/// </summary>
		public IterationResult CurrentSolution => this[^2];

		/// <summary>
		///     The results of the last solution (penultimate solved iteration [i - 2]).
		/// </summary>
		public IterationResult LastSolution => this[^3];

		#endregion

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
			: this(number, forces, Vector<double>.Build.Dense(forces.Count), Matrix<double>.Build.Dense(forces.Count, forces.Count), parameters)
		{
			LoadFactor = loadFactor;
		}

		/// <inheritdoc cref="LoadStep" />
		/// <param name="initialDisplacements">The initial displacement vector of this step.</param>
		/// <param name="stiffness">The stiffness matrix of this step.</param>
		/// <param name="parameters">The analysis parameters.</param>
		private LoadStep(int number, Vector<double> forces, Vector<double> initialDisplacements, Matrix<double> stiffness, AnalysisParameters parameters)
		{
			Number               = number;
			Forces               = forces;
			InitialDisplacements = initialDisplacements;
			Parameters           = parameters;
			
			Add(new IterationResult(initialDisplacements, Vector<double>.Build.Dense(initialDisplacements.Count), stiffness));
		}

		#endregion

		#region Methods

		/// <summary>
		///		Do the initial load step of a nonlinear analysis procedure.
		/// </summary>
		/// <param name="femInput">The finite element input.</param>
		/// <param name="loadFactor">The load factor of the first step. </param>
		/// <param name="parameters">The analysis parameters.</param>
		/// <returns>
		///		The initial <see cref="LoadStep"/>.
		/// </returns>
		public static LoadStep InitialStep(IFEMInput<IFiniteElement> femInput, double loadFactor, AnalysisParameters parameters)
		{
			var step      = new LoadStep(femInput.ForceVector * loadFactor, loadFactor, parameters, 1);
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
		
		/// <summary>
		///		Add a new iteration in this load step.
		/// </summary>
		///  <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public void NewIteration(bool simulate = false)
		{
			Add(
				this.Any()
					? this.Last().Clone()
					: IterationResult.FromStepResult(this, simulate));
			
			// Increase iteration count
			this.Last().Number++;
		}

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

			if (!monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = Length.FromMillimeters(FinalDisplacements[monitoredIndex.Value]);

			// Set to step
			MonitoredDisplacement = new MonitoredDisplacement(disp, LoadFactor);
		}

		/// <summary>
		///		Get the accumulated displacement increment from the beginning of this load step until a final index.
		/// </summary>
		/// <param name="finalIndex">The final index to consider increments.</param>
		public Vector<double> AccumulatedDisplacementIncrement(Index finalIndex) =>
			this[finalIndex].Displacements - InitialDisplacements;

		/// <summary>
		///		Get the total accumulated displacement increment at this load step.
		/// </summary>
		public Vector<double> AccumulatedDisplacementIncrement() => AccumulatedDisplacementIncrement(^1);

		/// <inheritdoc />
		public LoadStep Clone() => new(Number, Forces.Clone(), FinalDisplacements.Clone(), Stiffness.Clone(), Parameters)
		{
			LoadFactor = LoadFactor
		};

		/// <summary>
		///     Check if iterative procedure must stop by achieving convergence or achieving the maximum number of iterations.
		/// </summary>
		/// <remarks>
		///     If the maximum number of iterations is reached, <see cref="IterationResult.Stop" /> is set to true.
		/// </remarks>
		/// <returns>
		///     True if convergence is reached or the maximum number of iterations is reached.
		/// </returns>
		protected bool IterativeStop()
		{
			// Check if one stop condition is reached
			OngoingIteration.Stop = OngoingIteration >= Parameters.MaxIterations           || OngoingIteration.ResidualForces.ContainsNaNOrInfinity() ||
									OngoingIteration.Displacements.ContainsNaNOrInfinity() || OngoingIteration.Stiffness.ContainsNaN();

			return 
				OngoingIteration.Stop ||
				VerifyConvergence(OngoingIteration.ForceConvergence, Parameters.ForceTolerance) ||
				VerifyConvergence(OngoingIteration.DisplacementConvergence, Parameters.DisplacementTolerance);
		}
		
		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="convergence"> Calculated convergence. </param>
		/// <param name="tolerance">The required tolerance.</param>
		/// <param name="minIterations">The minimum number of iterations.</param>
		protected bool VerifyConvergence(double convergence, double tolerance, int minIterations = 2) => convergence <= tolerance && OngoingIteration >= minIterations;

		/// <inheritdoc />
		public override string ToString() => $"Load step {Number}";

		#endregion

		#endregion

		#region Operators

		/// <summary>
		///     Get the number of a step.
		/// </summary>
		/// <returns>
		///     <see cref="LoadStep.Number" />
		/// </returns>
		public static explicit operator int(LoadStep loadStep) => loadStep.Number;

		/// <summary>
		///		Check the step number.
		/// </summary>
		/// <returns>
		///		True if the step number is equal to the right number.
		/// </returns>
		public static bool operator ==(LoadStep left, int right) => left.Number == right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is not equal to the right number.
		/// </returns>
		public static bool operator !=(LoadStep left, int right) => left.Number != right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is bigger than the right number.
		/// </returns>
		public static bool operator >(LoadStep left, int right) => left.Number > right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is smaller than the right number.
		/// </returns>
		public static bool operator <(LoadStep left, int right) => left.Number < right;
		
		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is bigger or equal to the right number.
		/// </returns>
		public static bool operator >=(LoadStep left, int right) => left.Number >= right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is smaller or equal to the right number.
		/// </returns>
		public static bool operator <=(LoadStep left, int right) => left.Number <= right;

		#endregion

	}
}