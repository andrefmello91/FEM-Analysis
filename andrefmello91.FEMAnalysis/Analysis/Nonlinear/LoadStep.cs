using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.FEMAnalysis;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

using static andrefmello91.FEMAnalysis.Analysis;
using static andrefmello91.FEMAnalysis.StiffnessMatrix;
using static andrefmello91.FEMAnalysis.NonlinearAnalysis;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Generic class for load step.
	/// </summary>
	public class LoadStep: IEnumerable<IIteration>
	{
		/// <summary>
		///		The vector of full applied forces.
		/// </summary>
		protected readonly ForceVector FullForceVector;

		/// <summary>
		///		Auxiliary iteration list.
		/// </summary>
		protected readonly List<IIteration> Iterations = new();
		
		#region Properties

		/// <summary>
		///     The status of this step. True if convergence was reached.
		/// </summary>
		public bool Converged { get; protected set; }

		/// <summary>
		///     The results of the current solution (last solved iteration [i - 1]).
		/// </summary>
		public IIteration LastIteration => Iterations.Count > 1
			? Iterations[^2]
			: Iteration.From(FullForceVector.Count, this is SimulationStep);

		/// <summary>
		///     The displacement convergence of this step.
		/// </summary>
		public double DisplacementConvergence => CurrentIteration.DisplacementConvergence;

		/// <summary>
		///		The total displacement increment at this step.
		/// </summary>
		public DisplacementVector DisplacementIncrement => FinalDisplacements - InitialDisplacements;
		
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
		public ForceVector Forces => LoadFactor * FullForceVector;

		/// <summary>
		///     The displacement vector at the beginning of this step.
		/// </summary>
		public DisplacementVector InitialDisplacements { get; }

		/// <summary>
		///     The results of the last solution (penultimate solved iteration [i - 2]).
		/// </summary>
		public IIteration PenultimateIteration => Iterations.Count > 2
			? Iterations[^3]
			: Iteration.From(FullForceVector.Count, this is SimulationStep);

		/// <summary>
		///     The load factor of this step.
		/// </summary>
		public double LoadFactor { get; protected set; }

		/// <summary>
		///		The load factor increment of this step.
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
		///     The results of the ongoing iteration.
		/// </summary>
		public IIteration CurrentIteration => Iterations[^1];

		/// <summary>
		///     The analysis parameters.
		/// </summary>
		public AnalysisParameters Parameters { get; }

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

		/// <summary>
		///		Get the iteration at this index.
		/// </summary>
		public IIteration this[int index] => Iterations[index];
		
		/// <summary>
		///		Get the iteration at this index.
		/// </summary>
		public IIteration this[Index index] => Iterations[index];
		
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
			: this(number, fullForceVector, loadFactor, DisplacementVector.Zero(fullForceVector.Count), Zero(fullForceVector.Count), parameters, simulate)
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
			LoadFactor           = loadFactor;
			InitialDisplacements = initialDisplacements;
			Parameters           = parameters;
			
			Iterations.Add(Iteration.From(initialDisplacements, ForceVector.Zero(fullForceVector.Count), stiffness, simulate));
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
				false => new LoadStep(femInput.ForceVector, loadFactor, parameters, stepNumber, simulate),
				_     => new SimulationStep(femInput.ForceVector, loadFactor, parameters, stepNumber)
			};

		/// <summary>
		///     Do the initial load step of a nonlinear analysis procedure.
		/// </summary>
		/// <inheritdoc cref="From(ForceVector,double,DisplacementVector,StiffnessMatrix,AnalysisParameters,int,bool)"/>
		/// <returns>
		///     The initial <see cref="LoadStep" />.
		/// </returns>
		public static LoadStep InitialStep(IFEMInput femInput, AnalysisParameters parameters, bool simulate = false)
		{
			if (simulate)
				return SimulationStep.InitialStep(femInput, parameters);
			
			var lf   = StepIncrement(parameters.NumberOfSteps);
			
			var step = From(femInput, lf, parameters, 1, simulate);
			
			var iteration = step.CurrentIteration;

			// Get the initial stiffness and force vector simplified
			iteration.Stiffness = Assemble(femInput);

			// Calculate initial displacements
			iteration.IncrementDisplacements(iteration.Stiffness.Solve(step.Forces));

			// Update displacements in grips and elements
			femInput.Grips.SetDisplacements(iteration.Displacements);
			femInput.UpdateDisplacements();

			// Calculate element forces
			femInput.CalculateForces();

			// Update internal forces
			iteration.UpdateForces(step.Forces, ForceVector.AssembleInternal(femInput));

			return step;
		}

		///  <summary>
		/// 		Create a load step from the last load step.
		///  </summary>
		///  <param name="lastStep">The last load step.</param>
		///  <param name="incrementLoad">Increment load of the new step? If it's a <see cref="SimulationStep"/>, load is not increased.</param>
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
		///     Get the accumulated displacement increment from the beginning of this load step until a final index.
		/// </summary>
		/// <param name="finalIndex">The final index to consider increments.</param>
		public DisplacementVector AccumulatedDisplacementIncrement(Index finalIndex)
		{
			var iterations = Iterations.Where(i => i.Number > 0).ToList();
			
			return iterations.Count < finalIndex.Value
				?  DisplacementVector.Zero(InitialDisplacements.Count)
				: iterations[finalIndex].Displacements - InitialDisplacements;
		}

		/// <summary>
		///     Increment forces in this step by the default load factor increment.
		/// </summary>
		private void IncrementLoad() => IncrementLoad(StepIncrement(Parameters.NumberOfSteps));
		
		/// <summary>
		///     Increment forces in this step by a custom load factor increment.
		/// </summary>
		/// <param name="loadFactorIncrement">The increment of the load factor.</param>
		public void IncrementLoad(double loadFactorIncrement) => LoadFactor += loadFactorIncrement;

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
				// Add iteration
				NewIteration();

				// Update displacements, stiffness and forces
				UpdateDisplacements(femInput);
				UpdateStiffness();
				UpdateForces(femInput);
				
				// Calculate convergence
				CurrentIteration.CalculateConvergence(Forces, FirstIteration.DisplacementIncrement);
				
			} while (!IterativeStop());
		}

		/// <summary>
		///		Update displacements in elements and calculate internal forces.
		/// </summary>
		protected void UpdateElements(IFEMInput femInput)
		{
			// Update elements
			femInput.Grips.SetDisplacements(CurrentIteration.Displacements);
			femInput.UpdateDisplacements();
			femInput.CalculateForces();
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
		///		Update forces and calculate convergence.
		/// </summary>
		protected virtual void UpdateForces(IFEMInput femInput)
		{
			// Calculate element forces
			femInput.CalculateForces();

			// Update internal forces
			var intForces = ForceVector.AssembleInternal(femInput);
			CurrentIteration.UpdateForces(Forces, intForces);
		}
		
		/// <summary>
		///     Update displacements.
		/// </summary>
		private void UpdateDisplacements(IFEMInput femInput)
		{
			var curIt  = CurrentIteration;
			var lastIt = LastIteration;

			// Calculate increment from residual
			var dUr = -curIt.Stiffness.Solve(lastIt.ResidualForces);
			curIt.IncrementDisplacements(dUr);
			
			// Update displacements in grips and elements
			femInput.Grips.SetDisplacements(curIt.Displacements);
			femInput.UpdateDisplacements();
		}

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix{T}" /> of current iteration.
		/// </summary>
		protected void UpdateStiffness()
		{
			// If analysis stopped or solver is modified newton raphson and step didn't converge yet
			if (Stop || !Converged && Parameters.Solver is NonLinearSolver.ModifiedNewtonRaphson)
				return;

			CurrentIteration.Stiffness += StiffnessIncrement(CurrentIteration, LastIteration, Parameters.Solver);
		}


		/// <summary>
		///     Set step results after achieving convergence.
		/// </summary>
		public void SetResults(int? monitoredIndex = null)
		{
			if (!monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = FinalDisplacements[monitoredIndex.Value];

			// Set to step
			MonitoredDisplacement = new MonitoredDisplacement(disp, LoadFactor);
		}

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

		#region Interface Implementations

		#endregion

		#region Object override

		/// <inheritdoc />
		public IEnumerator<IIteration> GetEnumerator() => Iterations.GetEnumerator();

		/// <inheritdoc />
		public override string ToString() => $"Load step {Number}";

		/// <inheritdoc />
		IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

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
		///     Check the step number.
		/// </summary>
		/// <returns>
		///     True if the step number is equal to the right number.
		/// </returns>
		public static bool operator ==(LoadStep left, int right) => left.Number == right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is not equal to the right number.
		/// </returns>
		public static bool operator !=(LoadStep left, int right) => left.Number != right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is bigger than the right number.
		/// </returns>
		public static bool operator >(LoadStep left, int right) => left.Number > right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is smaller than the right number.
		/// </returns>
		public static bool operator <(LoadStep left, int right) => left.Number < right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is bigger or equal to the right number.
		/// </returns>
		public static bool operator >=(LoadStep left, int right) => left.Number >= right;

		/// <inheritdoc cref="op_Equality" />
		/// <returns>
		///     True if the step number is smaller or equal to the right number.
		/// </returns>
		public static bool operator <=(LoadStep left, int right) => left.Number <= right;

		#endregion

	}
}