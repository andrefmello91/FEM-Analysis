using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.OnPlaneComponents;
using MathNet.Numerics.LinearAlgebra;
using static andrefmello91.FEMAnalysis.NonlinearAnalysis;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Generic class for load step.
	/// </summary>
	public class LoadStep : IEnumerable<IIteration>
	{

		#region Fields

		/// <summary>
		///     The vector of full applied forces.
		/// </summary>
		/// <remarks>
		///     Simplified at constrained DoFs.
		/// </remarks>
		protected readonly ForceVector FullForceVector;

		/// <summary>
		///     Auxiliary iteration list.
		/// </summary>
		protected readonly List<IIteration> Iterations = new();

		private double _loadFactor;

		#endregion

		#region Properties

		/// <summary>
		///     Get the iteration at this index.
		/// </summary>
		public IIteration this[int index] => Iterations[index];

		/// <summary>
		///     Get the iteration at this index.
		/// </summary>
		public IIteration this[Index index] => Iterations[index];

		/// <summary>
		///     The status of this step. True if convergence was reached.
		/// </summary>
		public bool Converged { get; protected set; }

		/// <summary>
		///     The results of the ongoing iteration.
		/// </summary>
		public IIteration CurrentIteration => Iterations[^1];

		/// <summary>
		///     The displacement convergence of this step.
		/// </summary>
		public double DisplacementConvergence => CurrentIteration.DisplacementConvergence;

		/// <summary>
		///     The total displacement increment at this step.
		/// </summary>
		public DisplacementVector DisplacementIncrement => (DisplacementVector) (FinalDisplacements - InitialDisplacements);

		/// <summary>
		///     The displacement vector at the end of this step.
		/// </summary>
		public DisplacementVector FinalDisplacements => Iterations.Last().Displacements;

		/// <summary>
		///     Get the first iteration of the current step.
		/// </summary>
		public IIteration FirstIteration => Iterations.Find(i => i.Number == 1)!;

		/// <summary>
		///     The force convergence of this step.
		/// </summary>
		public double ForceConvergence => CurrentIteration.ForceConvergence;

		/// <summary>
		///     The force vector of this step.
		/// </summary>
		/// <inheritdoc cref="FullForceVector" />
		public ForceVector Forces => (ForceVector) (LoadFactor * FullForceVector);

		/// <summary>
		///     The displacement vector at the beginning of this step.
		/// </summary>
		public DisplacementVector InitialDisplacements { get; }

		/// <summary>
		///     The results of the current solution (last solved iteration [i - 1]).
		/// </summary>
		public IIteration LastIteration => Iterations.Count > 1
			? Iterations[^2]
			: Iteration.From(FullForceVector.Count, this is SimulationStep);

		/// <summary>
		///     The load factor of this step.
		/// </summary>
		public virtual double LoadFactor => _loadFactor;

		/// <summary>
		///     The load factor increment of this step.
		/// </summary>
		public virtual double LoadFactorIncrement => StepIncrement(Parameters.NumberOfSteps);

		/// <summary>
		///     The monitored displacement of this step.
		/// </summary>
		public MonitoredDisplacement? MonitoredDisplacement { get; protected set; }

		/// <summary>
		///     The number of this step.
		/// </summary>
		public int Number { get; set; }

		/// <summary>
		///     The analysis parameters.
		/// </summary>
		public AnalysisParameters Parameters { get; }

		/// <summary>
		///     The results of the last solution (penultimate solved iteration [i - 2]).
		/// </summary>
		public IIteration PenultimateIteration => Iterations.Count > 2
			? Iterations[^3]
			: Iteration.From(FullForceVector.Count, this is SimulationStep);

		/// <summary>
		///     The current stiffness matrix of this step (stiffness of the current iteration.
		/// </summary>
		public StiffnessMatrix Stiffness => Iterations.Last().Stiffness;

		/// <summary>
		///     Get/set when to stop analysis.
		/// </summary>
		/// <remarks>
		///     If true, convergence was not reached at this load step.
		/// </remarks>
		public bool Stop { get; private set; }

		#endregion

		#region Constructors

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="loadFactor">The load factor of this step.</param>
		/// <param name="number">The number of this step.</param>
		/// <param name="fullForceVector">The full applied force vector of the model.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		protected LoadStep(ForceVector fullForceVector, double loadFactor, AnalysisParameters parameters, int number = 0, bool simulate = false)
			: this(number, fullForceVector, loadFactor, DisplacementVector.Zero(fullForceVector.Count), StiffnessMatrix.Zero(fullForceVector.Count), parameters, simulate)
		{
		}

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="initialDisplacements">The initial displacement vector of this step.</param>
		/// <param name="stiffness">The stiffness matrix of this step.</param>
		/// <param name="parameters">The analysis parameters.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		protected LoadStep(int number, ForceVector fullForceVector, double loadFactor, DisplacementVector initialDisplacements, StiffnessMatrix stiffness, AnalysisParameters parameters, bool simulate = false)
		{
			Number               = number;
			FullForceVector      = fullForceVector;
			_loadFactor          = loadFactor;
			InitialDisplacements = initialDisplacements;
			Parameters           = parameters;

			Iterations.Add(Iteration.From(initialDisplacements, ForceVector.Zero(fullForceVector.Count), stiffness, simulate, loadFactor));
		}

		#endregion

		#region Methods

		/// <summary>
		///     Create a load step for a nonlinear analysis procedure.
		/// </summary>
		/// <param name="fullForceVector">The full applied force vector of the model.</param>
		/// <param name="loadFactor">The load factor of this step.</param>
		/// <param name="initialDisplacements">The initial displacement vector of this step.</param>
		/// <param name="stiffness">The stiffness matrix of this step.</param>
		/// <param name="parameters">The analysis parameters.</param>
		/// <param name="stepNumber">The number of the load step.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public static LoadStep From(ForceVector fullForceVector, double loadFactor, DisplacementVector initialDisplacements, StiffnessMatrix stiffness, AnalysisParameters parameters, int stepNumber, bool simulate = false) =>
			simulate switch
			{
				false => new LoadStep(stepNumber, fullForceVector, loadFactor, initialDisplacements, stiffness, parameters, simulate),
				_     => new SimulationStep(stepNumber, fullForceVector, loadFactor, initialDisplacements, stiffness, parameters)
			};



		/// <summary>
		///     Create a load step for a nonlinear analysis procedure.
		/// </summary>
		/// <param name="femInput">The finite element input.</param>
		/// <param name="loadFactor">The load factor of the first step. </param>
		/// <param name="parameters">The analysis parameters.</param>
		/// <param name="stepNumber">The number of the load step.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public static LoadStep From(IFEMInput femInput, double loadFactor, AnalysisParameters parameters, int stepNumber, bool simulate = false) =>
			simulate switch
			{
				false => new LoadStep(femInput.AssembleExternalForces(), loadFactor, parameters, stepNumber, simulate),
				_     => new SimulationStep(femInput.AssembleExternalForces(), loadFactor, parameters, stepNumber)
			};

		/// <summary>
		///     Create a load step from the last load step.
		/// </summary>
		/// <param name="lastStep">The last load step.</param>
		/// <param name="incrementLoad">
		///     Increment load of the new step? If it's a <see cref="SimulationStep" />, load is not
		///     increased.
		/// </param>
		public static LoadStep FromLastStep(LoadStep lastStep, bool incrementLoad = true)
		{
			if (lastStep is SimulationStep simulationStep)
				return SimulationStep.FromLastStep(simulationStep);

			var newStep = From(lastStep.FullForceVector, lastStep.LoadFactor, lastStep.FinalDisplacements, lastStep.Stiffness, lastStep.Parameters, lastStep.Number + 1);

			if (incrementLoad)
				newStep.IncrementLoad();

			return newStep;
		}

		/// <summary>
		///     Do the initial load step of a nonlinear analysis procedure.
		/// </summary>
		/// <inheritdoc cref="From(ForceVector,double,DisplacementVector,StiffnessMatrix,AnalysisParameters,int,bool)" />
		/// <returns>
		///     The initial <see cref="LoadStep" />.
		/// </returns>
		public static LoadStep InitialStep(IFEMInput femInput, AnalysisParameters parameters)
		{
			var lf = StepIncrement(parameters.NumberOfSteps);

			var step = From(femInput, lf, parameters, 1);

			var iteration = step.CurrentIteration;

			// Get the initial stiffness and force vector simplified
			iteration.Stiffness = femInput.AssembleStiffness();

			// Update internal forces
			iteration.UpdateForces(step.Forces, ForceVector.Zero(femInput.NumberOfDoFs));

			return step;
		}

		/// <summary>
		///     Get the accumulated displacement increment from the beginning of this load step until a final index.
		/// </summary>
		/// <param name="finalIndex">The final index to consider increments.</param>
		public DisplacementVector AccumulatedDisplacementIncrement(Index finalIndex)
		{
			var iterations = Iterations.Where(i => i.Number > 0).ToArray();

			DisplacementVector accD;

			try
			{
				accD = (DisplacementVector) (iterations[finalIndex].Displacements - InitialDisplacements);
			}
			catch
			{
				accD = DisplacementVector.Zero(InitialDisplacements.Count);
			}

			return accD;
		}

		/// <summary>
		///     Increment forces in this step by a custom load factor increment.
		/// </summary>
		/// <param name="loadFactorIncrement">The increment of the load factor.</param>
		public void IncrementLoad(double loadFactorIncrement) => _loadFactor += loadFactorIncrement;

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		/// <param name="femInput">The finite element input.</param>
		public virtual void Iterate(IFEMInput femInput)
		{
			// Initiate first iteration
			foreach (var iteration in Iterations)
				iteration.Number = 0;

			// Iterate
			do
			{
				// Update stiffness
				UpdateStiffness();

				// Add iteration
				NewIteration();

				// Update displacements
				UpdateDisplacements();

				// Update elements and forces
				UpdateElements(femInput);
				UpdateForces(femInput);

				// Calculate convergence
				CurrentIteration.CalculateConvergence(Forces, FirstIteration.DisplacementIncrement);

			} while (!IterativeStop());

			if (!Stop && Parameters.Solver is NonLinearSolver.ModifiedNewtonRaphson)
				UpdateStiffness();
		}


		/// <summary>
		///     Set step results after achieving convergence.
		/// </summary>
		public void SetResults(IFEMInput femInput, int? monitoredIndex = null)
		{
			if (!monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = FinalDisplacements[monitoredIndex.Value];

			// Set to step
			MonitoredDisplacement = new MonitoredDisplacement(disp, LoadFactor);

			foreach (var element in femInput.MonitoredElements)
				element.AddValue(LoadFactor);
		}

		/// <inheritdoc />
		public override string ToString() => $"Load step {Number}";

		/// <summary>
		///     Check if iterative procedure must stop by achieving convergence or achieving the maximum number of iterations.
		/// </summary>
		/// <returns>
		///     True if convergence is reached or the maximum number of iterations is reached.
		/// </returns>
		protected bool IterativeStop()
		{
			// Check if one stop condition is reached
			Stop = CurrentIteration.CheckStopCondition(Parameters);

			// Check convergence
			Converged = CurrentIteration.CheckConvergence(Parameters);

			return
				Stop || Converged;
		}

		/// <summary>
		///     Add a new iteration in this load step.
		/// </summary>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		protected void NewIteration(bool simulate = false)
		{
			if (CurrentIteration.Number > 0)
				Iterations.Add(this.Any()
					? this.Last().Clone()
					: Iteration.FromStepResult(this, simulate));

			// Increase iteration count
			CurrentIteration.Number++;
		}

		/// <summary>
		///     Update displacements in elements and calculate internal forces.
		/// </summary>
		protected void UpdateElements(IFEMInput femInput)
		{
			// Update elements
			femInput.Grips.SetDisplacements(CurrentIteration.Displacements);
			femInput.UpdateDisplacements();
			femInput.CalculateForces();
		}


		/// <summary>
		///     Update forces and calculate convergence.
		/// </summary>
		protected void UpdateForces(IFEMInput femInput) => CurrentIteration.UpdateForces(Forces, femInput.AssembleInternalForces());

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix{T}" /> of current iteration.
		/// </summary>
		protected void UpdateStiffness()
		{
			// If analysis stopped or solver is modified newton raphson and step didn't converge yet
			if (CurrentIteration.Number <= 1 || Stop || !Converged && Parameters.Solver is NonLinearSolver.ModifiedNewtonRaphson)
				return;

			CurrentIteration.Stiffness = (StiffnessMatrix) (CurrentIteration.Stiffness + StiffnessIncrement(CurrentIteration, LastIteration, Parameters.Solver));
		}

		/// <summary>
		///     Increment forces in this step by the default load factor increment.
		/// </summary>
		private void IncrementLoad() => IncrementLoad(StepIncrement(Parameters.NumberOfSteps));

		/// <summary>
		///     Update displacements.
		/// </summary>
		private void UpdateDisplacements()
		{
			var curIt  = CurrentIteration;
			var lastIt = LastIteration;

			// Calculate increment from residual
			var dUr = (DisplacementVector) (-curIt.Stiffness.Solve(lastIt.ResidualForces));
			curIt.IncrementDisplacements(dUr);

			// Update displacements in grips and elements
			// femInput.Grips.SetDisplacements(curIt.Displacements);
			// femInput.UpdateDisplacements();
		}

		/// <inheritdoc />
		public IEnumerator<IIteration> GetEnumerator() => Iterations.GetEnumerator();

		/// <inheritdoc />
		IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

		#endregion

		#region Operators

		/// <summary>
		///     Check the step number.
		/// </summary>
		/// <returns>
		///     True if the step number is equal to the right number.
		/// </returns>
		public static bool operator ==(LoadStep left, int right) => left.Number == right;

		/// <summary>
		///     Get the number of a step.
		/// </summary>
		/// <returns>
		///     <see cref="LoadStep.Number" />
		/// </returns>
		public static explicit operator int(LoadStep loadStep) => loadStep.Number;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is bigger than the right number.
		/// </returns>
		public static bool operator >(LoadStep left, int right) => left.Number > right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is bigger or equal to the right number.
		/// </returns>
		public static bool operator >=(LoadStep left, int right) => left.Number >= right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is not equal to the right number.
		/// </returns>
		public static bool operator !=(LoadStep left, int right) => left.Number != right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is smaller than the right number.
		/// </returns>
		public static bool operator <(LoadStep left, int right) => left.Number < right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is smaller or equal to the right number.
		/// </returns>
		public static bool operator <=(LoadStep left, int right) => left.Number <= right;

		#endregion

	}
}