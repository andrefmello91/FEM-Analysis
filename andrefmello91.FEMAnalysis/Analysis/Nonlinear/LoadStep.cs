using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.FEMAnalysis;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

using static andrefmello91.FEMAnalysis.Analysis<andrefmello91.FEMAnalysis.IFiniteElement>;
using static andrefmello91.FEMAnalysis.NonlinearAnalysis;
using static andrefmello91.FEMAnalysis.Simulation;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Generic class for load step.
	/// </summary>
	public class LoadStep: IEnumerable<IIteration>
	{
		/// <summary>
		///		True if this is a simulation.
		/// </summary>
		private bool _simulate;

		/// <summary>
		///		The vector of full applied forces.
		/// </summary>
		private readonly Vector<double> _fullForceVector;

		/// <summary>
		///		Auxiliary iteration list.
		/// </summary>
		private readonly List<IIteration> _iterations = new();
		
		#region Properties

		/// <summary>
		///     The status of this step. True if convergence was reached.
		/// </summary>
		public bool Converged { get; private set; }

		/// <summary>
		///     The results of the current solution (last solved iteration [i - 1]).
		/// </summary>
		public IIteration CurrentSolution => _iterations[^2];

		/// <summary>
		///     The displacement convergence of this step.
		/// </summary>
		public double DisplacementConvergence => OngoingIteration.DisplacementConvergence;

		/// <summary>
		///     The displacement vector at the end of this step.
		/// </summary>
		public Vector<double> FinalDisplacements => _iterations.Last().Displacements;

		/// <summary>
		///     Get the first iteration of the current step.
		/// </summary>
		public IIteration FirstIteration => _iterations.Find(i => i.Number == 1)!;

		/// <summary>
		///     The force convergence of this step.
		/// </summary>
		public double ForceConvergence => OngoingIteration.ForceConvergence;

		/// <summary>
		///     The force vector of this step.
		/// </summary>
		public Vector<double> Forces => LoadFactor * _fullForceVector;

		/// <summary>
		///     The displacement vector at the beginning of this step.
		/// </summary>
		public Vector<double> InitialDisplacements { get; }

		/// <summary>
		///     The results of the last solution (penultimate solved iteration [i - 2]).
		/// </summary>
		public IIteration LastSolution => _iterations[^3];

		/// <summary>
		///     The load factor of this step.
		/// </summary>
		public double LoadFactor { get; private set; }

		/// <summary>
		///     The monitored displacement of this step.
		/// </summary>
		public MonitoredDisplacement? MonitoredDisplacement { get; private set; }

		/// <summary>
		///     The number of this step.
		/// </summary>
		public int Number { get; set; }

		/// <summary>
		///     The results of the ongoing iteration.
		/// </summary>
		public IIteration OngoingIteration => _iterations[^1];

		/// <summary>
		///     The analysis parameters.
		/// </summary>
		public AnalysisParameters Parameters { get; }

		/// <summary>
		///     The current stiffness matrix of this step (stiffness of the current iteration.
		/// </summary>
		public Matrix<double> Stiffness => _iterations.Last().Stiffness;

		/// <summary>
		///		The required number of iterations for achieving convergence for this load step.
		/// </summary>
		public int RequiredIterations => _iterations.Count(i => i.Number > 0);

		/// <summary>
		///		The desired number of iterations for achieving convergence for this load step.
		/// </summary>
		/// <remarks>
		///		Default : 5
		///	</remarks>
		public int DesiredIterations { get; set; } = 5;
		
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
		public IIteration this[int index] => _iterations[index];
		
		/// <summary>
		///		Get the iteration at this index.
		/// </summary>
		public IIteration this[Index index] => _iterations[index];
		
		#endregion

		#region Constructors

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="loadFactor">The load factor of this step.</param>
		/// <param name="number">The number of this step.</param>
		/// <param name="fullForceVector">The force vector of this step.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		private LoadStep(Vector<double> fullForceVector, double loadFactor, AnalysisParameters parameters, int number = 0, bool simulate = false)
			: this(number, fullForceVector, loadFactor, Vector<double>.Build.Dense(fullForceVector.Count), Matrix<double>.Build.Dense(fullForceVector.Count, fullForceVector.Count), parameters, simulate)
		{
		}

		/// <inheritdoc cref="LoadStep" />
		/// <param name="initialDisplacements">The initial displacement vector of this step.</param>
		/// <param name="stiffness">The stiffness matrix of this step.</param>
		/// <param name="parameters">The analysis parameters.</param>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		private LoadStep(int number, Vector<double> fullForceVector, double loadFactor, Vector<double> initialDisplacements, Matrix<double> stiffness, AnalysisParameters parameters, bool simulate = false)
		{
			_simulate            = simulate;
			Number               = number;
			_fullForceVector     = fullForceVector;
			LoadFactor           = loadFactor;
			InitialDisplacements = initialDisplacements;
			Parameters           = parameters;
			
			_iterations.Add(Iteration.From(initialDisplacements, Vector<double>.Build.Dense(initialDisplacements.Count), stiffness, simulate));
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
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		public static LoadStep From(IFEMInput<IFiniteElement> femInput, double loadFactor, AnalysisParameters parameters, int stepNumber, bool simulate = false) =>
			new(femInput.ForceVector, loadFactor, parameters, stepNumber, simulate);

		/// <summary>
		///     Do the initial load step of a nonlinear analysis procedure.
		/// </summary>
		/// <inheritdoc cref="From"/>
		/// <returns>
		///     The initial <see cref="LoadStep" />.
		/// </returns>
		public static LoadStep InitialStep(IFEMInput<IFiniteElement> femInput, AnalysisParameters parameters, bool simulate = false)
		{
			var lf   = StepIncrement(parameters.NumberOfSteps);
			
			var step = From(femInput, lf, parameters, 1, simulate);
			
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
			iteration.UpdateForces(fi, femInput.AssembleInternalForces());

			return step;
		}

		///  <summary>
		/// 		Create a load step from the last load step.
		///  </summary>
		///  <param name="lastStep">The last load step.</param>
		///  <param name="incrementLoad">Increment load of the new step?</param>
		///  <remarks>
		/// 		This method doesn't increase load, only the step number.
		///  </remarks>
		public static LoadStep FromLastStep(LoadStep lastStep, bool incrementLoad = true)
		{
			var newStep = new LoadStep(lastStep.Number + 1, lastStep._fullForceVector, lastStep.LoadFactor, lastStep.FinalDisplacements, lastStep.Stiffness, lastStep.Parameters, lastStep._simulate);
			
			if (incrementLoad)
				newStep.IncrementLoad();

			return newStep;
		}

		/// <summary>
		///     Get the accumulated displacement increment from the beginning of this load step until a final index.
		/// </summary>
		/// <param name="finalIndex">The final index to consider increments.</param>
		public Vector<double> AccumulatedDisplacementIncrement(Index finalIndex) =>
			_iterations[finalIndex].Displacements - InitialDisplacements;

		/// <summary>
		///     Get the total accumulated displacement increment at this load step.
		/// </summary>
		public Vector<double> AccumulatedDisplacementIncrement() => AccumulatedDisplacementIncrement(^1);

		/// <summary>
		///     Increment forces in this step by the default load factor increment.
		/// </summary>
		public void IncrementLoad() => IncrementLoad(StepIncrement(Parameters.NumberOfSteps));
		
		/// <summary>
		///     Increment forces in this step by a custom load factor increment.
		/// </summary>
		/// <param name="loadFactorIncrement">The increment of the load factor.</param>
		public void IncrementLoad(double loadFactorIncrement)
		{
			if (_simulate && _iterations.Any() && _iterations.Last() is SimulationIteration itResult)
				itResult.LoadFactorIncrement = loadFactorIncrement;

			// Update values
			LoadFactor += loadFactorIncrement;
		}

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		/// <param name="femInput">The finite element input.</param>
		public virtual void Iterate(IFEMInput<IFiniteElement> femInput)
		{
			// Initiate first iteration
			foreach (var iteration in _iterations)
				iteration.Number = 0;

			// Iterate
			do
			{
				// Add iteration
				NewIteration(_simulate);

				// Do the initial iteration
				if(_simulate && OngoingIteration.Number == 1)
					InitialIteration(femInput);
				
				else
				{
					// Update stiffness and displacements
					UpdateDisplacements(femInput);
					UpdateStiffness(femInput);
				}

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

		/// <summary>
		///     Add a new iteration in this load step.
		/// </summary>
		/// <param name="simulate">Set true if the performed analysis is a simulation.</param>
		private void NewIteration(bool simulate = false)
		{
			_iterations.Add(this.Any()
				? this.Last().Clone()
				: Iteration.FromStepResult(this, simulate));

			// Increase iteration count
			OngoingIteration.Number++;
		}

		/// <summary>
		///		Steps to perform at the initial iteration of a simulation.
		/// </summary>
		private void InitialIteration(IFEMInput<IFiniteElement> femInput)
		{
			if (!_simulate)
				return;
			
			var stiffness = SimplifiedStiffness(OngoingIteration.Stiffness, femInput.ConstraintIndex);
			var ongIt     = (SimulationIteration) OngoingIteration;
			
			switch (Number)
			{
				// First iteration of first load step
				case 1:
					// Set initial increment
					ongIt.LoadFactorIncrement = StepIncrement(Parameters.NumberOfSteps);

					// Set initial residual
					var intForces = stiffness * ongIt.Displacements;
					ongIt.UpdateForces(_fullForceVector, intForces);

					// Calculate the initial increments
					var dUr = -stiffness.Solve(ongIt.ResidualForces);
					var dUf = stiffness.Solve(_fullForceVector);
					ongIt.IncrementDisplacements(dUr, dUf);
					
					// Calculate arc length
					CalculateArcLength(this);
					
					break;
				
				// First iteration of any load step except the first
				default:
					// Check increment sign
					ongIt.CalculateStiffnessParameter(this);
					
					// Calculate increments
					var rInc = -stiffness.Solve(CurrentSolution.ResidualForces);
					var fInc =  stiffness.Solve(SimplifiedForces(_fullForceVector, femInput.ConstraintIndex));
					
					// Set increments
					ongIt.LoadFactorIncrement = StepIncrement(this);
					ongIt.IncrementDisplacements(rInc, fInc);
					
					break;
			}
			
			// Update displacements in grips and elements
			femInput.Grips.SetDisplacements(ongIt.Displacements);
			femInput.UpdateDisplacements();
		}
		
		/// <summary>
		///     Update displacements.
		/// </summary>
		private void UpdateDisplacements(IFEMInput<IFiniteElement> femInput)
		{
			var ongIt  = OngoingIteration;
			var curSol = CurrentSolution;

			// Increment displacements
			var stiffness = SimplifiedStiffness(ongIt.Stiffness, femInput.ConstraintIndex);
			
			// Calculate increment from residual
			var dUr = -stiffness.Solve(curSol.ResidualForces);
			
			if (!_simulate)
				ongIt.IncrementDisplacements(dUr);
			
			else
			{
				var simIt = (SimulationIteration) ongIt;
				
				// Calculate increment from external forces
				var dUf = stiffness.Solve(_fullForceVector);
				
				simIt.IncrementDisplacements(dUr, dUf);
			}
				
			// Update displacements in grips and elements
			femInput.Grips.SetDisplacements(ongIt.Displacements);
			femInput.UpdateDisplacements();
		}

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix{T}" /> of current iteration.
		/// </summary>
		private void UpdateStiffness(IFEMInput<IFiniteElement> femInput)
		{
			var ongIt = OngoingIteration;

			switch (Parameters.Solver)
			{
				case NonLinearSolver.Secant:
					// Increment current stiffness
					var curSol  = CurrentSolution;
					var lastSol = LastSolution;

					ongIt.Stiffness += SecantIncrement(curSol.Stiffness, curSol.Displacements, lastSol.Displacements, curSol.ResidualForces, lastSol.ResidualForces);
					break;

				// For Newton-Raphson
				case NonLinearSolver.NewtonRaphson:
				case NonLinearSolver.ModifiedNewtonRaphson when ongIt.Number == 1:
					// Update stiffness in elements
					femInput.UpdateStiffness();

					// Set new values
					ongIt.Stiffness = femInput.AssembleStiffness();

					break;

				default:
					return;
			}
		}


		/// <summary>
		///     Set step results after achieving convergence.
		/// </summary>
		public void SetResults(int? monitoredIndex = null)
		{
			if (!monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = Length.FromMillimeters(FinalDisplacements[monitoredIndex.Value]);

			// Set to step
			MonitoredDisplacement = new MonitoredDisplacement(disp, LoadFactor);
		}

		/// <summary>
		///     Check if iterative procedure must stop by achieving convergence or achieving the maximum number of iterations.
		/// </summary>
		/// <returns>
		///     True if convergence is reached or the maximum number of iterations is reached.
		/// </returns>
		private bool IterativeStop()
		{
			// Check if one stop condition is reached
			Stop = OngoingIteration.CheckStopCondition(Parameters);

			// Check convergence
			Converged = OngoingIteration.CheckConvergence(Parameters);

			return
				Stop || Converged;
		}

		#region Interface Implementations

		#endregion

		#region Object override

		/// <inheritdoc />
		public IEnumerator<IIteration> GetEnumerator() => _iterations.GetEnumerator();

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