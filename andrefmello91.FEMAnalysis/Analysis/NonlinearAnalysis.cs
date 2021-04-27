using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Nonlinear analysis class
	/// </summary>
	public class NonlinearAnalysis : Analysis<IFiniteElement>
	{

		#region Fields

		/// <summary>
		///     The list of load step results.
		/// </summary>
		private readonly List<LoadStepResult> _loadSteps = new();
		
		/// <summary>
		///     The list of iteration results.
		/// </summary>
		/// <remarks>
		///		This is cleaned at the beginning of a load step.
		/// </remarks>
		private readonly List<IterationResult> _iterations = new();

		/// <summary>
		///     The results of the current (ongoing) iteration.
		/// </summary>
		private IterationResult CurrentIteration => _iterations[^1];

		/// <summary>
		///     The results of the last iteration [i - 1]).
		/// </summary>
		private IterationResult LastIteration => _iterations[^2];
		
		/// <summary>
		///     The results of the penultimate iteration [i - 2].
		/// </summary>
		private IterationResult PenultimateIteration => _iterations[^3];

		/// <summary>
		///     Field to store the DoF index for monitored displacements.
		/// </summary>
		private int? _monitoredIndex;

		#endregion

		#region Properties

		/// <inheritdoc />
		/// <remarks>
		///     The displacements of current iteration.
		/// </remarks>
		public override Vector<double>? DisplacementVector
		{
			get => CurrentIteration.Displacements;
			protected set
			{
				if (value is null)
					return;

				// Update last iteration
				// LastIteration.Displacements = CurrentIteration.Displacements;

				CurrentIteration.Displacements = value;
			}
		}

		/// <inheritdoc />
		/// <remarks>
		///     The stiffness of current iteration.
		/// </remarks>
		public override Matrix<double>? GlobalStiffness
		{
			get => CurrentIteration.Stiffness;
			protected set
			{
				if (value is null)
					return;

				// Update last iteration
				// LastIteration.Stiffness = CurrentIteration.Stiffness;

				CurrentIteration.Stiffness = value;
			}
		}

		/// <summary>
		///     Get/set the maximum number of iterations.
		/// </summary>
		public int MaxIterations { get; set; }

		/// <summary>
		///     Get/set the minimum number of iterations.
		/// </summary>
		public int MinIterations { get; set; }

		/// <summary>
		///     Get/set the number of load steps to execute.
		/// </summary>
		public int NumLoadSteps { get; set; }

		/// <summary>
		///     The nonlinear equation solver.
		/// </summary>
		public NonLinearSolver Solver { get; set; }

		/// <summary>
		///     Get/set when to stop analysis.
		/// </summary>
		public bool Stop { get; private set; }

		/// <summary>
		///     Get/set the stop message.
		/// </summary>
		public string StopMessage { get; private set; } = string.Empty;

		/// <summary>
		///     Get/set the convergence tolerance.
		/// </summary>
		public double Tolerance { get; set; }

		/// <summary>
		///     The current load step result.
		/// </summary>
		private LoadStepResult CurrentLoadStep => _loadSteps.Last();

		/// <summary>
		///     The last load step result.
		/// </summary>
		private LoadStepResult LastLoadStep => _loadSteps.Count > 1
			? _loadSteps[^2]
			: CurrentLoadStep;

		/// <summary>
		///     Get/set the residual force vector of current iteration.
		/// </summary>
		private Vector<double> ResidualForces
		{
			get => CurrentIteration.ResidualForces;
			set
			{
				// Update last iteration
				// LastIteration.ResidualForces = CurrentIteration.ResidualForces;

				CurrentIteration.ResidualForces = value;
			}
		}

		#endregion

		#region Constructors

		/// <summary>
		///     Nonlinear analysis constructor.
		/// </summary>
		/// <param name="nonlinearInput">The <see cref="IFEMInput{INonlinearElement}" />.</param>
		/// <param name="solver">The <see cref="NonLinearSolver" /> to use.</param>
		/// <param name="numLoadSteps">The number of load steps to perform (default: 50).</param>
		/// <param name="tolerance">The convergence tolerance (default: 1E-6).</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step (default: 10000).</param>
		/// <param name="minIterations">Minimum number of iterations for each load step (default: 2).</param>
		public NonlinearAnalysis(
			IFEMInput<IFiniteElement> nonlinearInput,
			NonLinearSolver solver = NonLinearSolver.NewtonRaphson,
			int numLoadSteps = 50,
			double tolerance = 1E-6,
			int maxIterations = 10000,
			int minIterations = 2)
			: base(nonlinearInput)
		{
			Solver        = solver;
			NumLoadSteps  = numLoadSteps;
			Tolerance     = tolerance;
			MaxIterations = maxIterations;
			MinIterations = minIterations;
		}

		#endregion

		#region Methods

		/// <summary>
		///     Calculate the secant stiffness increment.
		/// </summary>
		/// <param name="currentStiffness">The stiffness matrix from current iteration.</param>
		/// <param name="currentDisplacements">The displacement vector from current iteration.</param>
		/// <param name="lastDisplacements">The displacement vector from the last iteration.</param>
		/// <param name="currentResidual">The residual force vector from current iteration.</param>
		/// <param name="lastResidual">The residual force vector from last iteration.</param>
		/// <returns>
		///     <see cref="Matrix{T}" />
		/// </returns>
		public static Matrix<double> SecantIncrement(Matrix<double> currentStiffness, Vector<double> currentDisplacements, Vector<double> lastDisplacements, Vector<double> currentResidual, Vector<double> lastResidual)
		{
			// Calculate the variation of displacements and residual as vectors
			Vector<double>
				dDisp = currentDisplacements - lastDisplacements,
				dRes  = currentResidual - lastResidual;

			return
				((dRes - currentStiffness * dDisp) / dDisp.Norm(2)).ToColumnMatrix() * dDisp.ToRowMatrix();
		}

		/// <summary>
		///     Calculate the tangent stiffness increment.
		/// </summary>
		/// <param name="currentStiffness">The stiffness matrix from current iteration.</param>
		/// <param name="lastStiffness">The stiffness matrix from last iteration.</param>
		/// <param name="currentDisplacements">The displacement vector from current iteration.</param>
		/// <param name="lastDisplacements">The displacement vector from the last iteration.</param>
		/// <returns>
		///     <see cref="Matrix{T}" />
		/// </returns>
		public static Matrix<double> TangentIncrement(Matrix<double> currentStiffness, Matrix<double> lastStiffness, Vector<double> currentDisplacements, Vector<double> lastDisplacements)
		{
			// Get displacement variation
			var du = currentDisplacements - lastDisplacements;

			// Get stiffness variation
			var dk = currentStiffness - lastStiffness;

			var inc = Matrix<double>.Build.Dense(dk.RowCount, dk.ColumnCount);

			// Increment elements of stiffness matrix
			for (var i = 0; i < inc.RowCount; i++)
			for (var j = 0; j < inc.ColumnCount; j++)
				if (du[j] != 0)
					for (var k = 0; k < inc.ColumnCount; k++)
						inc[i, j] += dk[i, k] / du[j] * currentDisplacements[k];

			return inc;
		}

		/// <summary>
		///     Calculate the force based convergence.
		/// </summary>
		/// <param name="residualForces">The residual forces of the current iteration.</param>
		/// <param name="appliedForces">The applied forces of the current load step.</param>
		private static double Convergence(IEnumerable<double> residualForces, IEnumerable<double> appliedForces)
		{
			double
				num = residualForces.Sum(n => n * n),
				den = 1 + appliedForces.Sum(n => n * n);

			return
				num / den;
		}

		/// <summary>
		///     Execute the analysis.
		/// </summary>
		/// <inheritdoc cref="Initiate" />
		/// <inheritdoc cref="StepAnalysis" />
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis{TFiniteElement}.ForceVector" /> (default: 1).</param>
		public void Execute(int? monitoredIndex = null, bool simulate = false, double loadFactor = 1)
		{
			// Get force vector
			ForceVector = FemInput.ForceVector * loadFactor;

			// Initiate lists
			Initiate(monitoredIndex);

			// Analysis by load steps
			StepAnalysis(simulate);

			// Set Reactions
			FemInput.Grips.SetReactions(GetReactions());
		}

		/// <summary>
		///     Generate an <see cref="FEMOutput" /> from analysis results.
		/// </summary>
		/// <returns>
		///     null if no monitored index was provided.
		/// </returns>
		public FEMOutput? GenerateOutput() =>
			!_monitoredIndex.HasValue
				? null
				: new FEMOutput(_loadSteps
					.Where(ls => ls.IsCalculated)
					.Select(ls => ls.MonitoredDisplacement!.Value)
					.ToList());

		/// <summary>
		///     Calculate the secant stiffness <see cref="Matrix{T}" /> of current iteration.
		/// </summary>
		/// <param name="simplify">Simplify stiffness?</param>
		protected override void UpdateStiffness(bool simplify = true)
		{
			switch (Solver)
			{
				case NonLinearSolver.Secant:
					// Increment current stiffness
					GlobalStiffness = LastIteration.Stiffness + SecantIncrement(LastIteration.Stiffness, LastIteration.Displacements, PenultimateIteration.Displacements, LastIteration.ResidualForces, PenultimateIteration.ResidualForces);
					break;

				// For Newton-Raphson
				default:
					// Update stiffness in elements
					FemInput.Elements.UpdateStiffness();

					// Set new values
					GlobalStiffness = FemInput.AssembleStiffness();

					break;
			}

			// Simplify
			if (simplify)
				Simplify(GlobalStiffness, null, FemInput.ConstraintIndex);
		}

		/// <summary>
		///     Correct results from last load step after not achieving convergence.
		/// </summary>
		private void CorrectResults()
		{
			// Set displacements from last load step
			DisplacementVector = LastLoadStep.Displacements;
			GlobalStiffness    = LastLoadStep.Stiffness;
			FemInput.Grips.SetDisplacements(DisplacementVector);
			FemInput.Elements.UpdateDisplacements();

			// Calculate element forces
			FemInput.Elements.CalculateForces();
		}

		/// <summary>
		///     Update displacements.
		/// </summary>
		private void DisplacementUpdate()
		{
			// Increment displacements
			DisplacementVector = LastIteration.Displacements - GlobalStiffness!.Solve(LastIteration.ResidualForces);

			// Update displacements in grips and elements
			FemInput.Grips.SetDisplacements(DisplacementVector);
			FemInput.Elements.UpdateDisplacements();
		}

		/// <summary>
		///     Initiate fields.
		/// </summary>
		/// <param name="monitoredIndex">The index of a degree of freedom to monitor, if wanted.</param>
		private void Initiate(int? monitoredIndex)
		{
			_monitoredIndex = monitoredIndex;

			// Initiate solution values
			for (int i = 0; i < 3 ; i++)
				_iterations.Add(new IterationResult(FemInput.NumberOfDoFs));

			// Get the initial stiffness and force vector simplified
			GlobalStiffness = FemInput.AssembleStiffness();
			Simplify(GlobalStiffness, ForceVector, FemInput.ConstraintIndex);

			// Calculate initial displacements
			var fi = ForceVector / NumLoadSteps;
			DisplacementVector = GlobalStiffness!.Solve(fi);

			// Add iterations
			// for (int i = 0; i < 2 ; i++)
			// 	_iterations.Add(CurrentIteration.Clone());
			
			// Update displacements in grips and elements
			FemInput.Grips.SetDisplacements(DisplacementVector);
			FemInput.Elements.UpdateDisplacements();
		}

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		private void Iterate()
		{
			// Clear iteration list
			if ((int) CurrentLoadStep > 1)
			{
				var lastIt = LastIteration.Clone();
				var curIt  = CurrentIteration.Clone();
				_iterations.Clear();
				_iterations.AddRange(new []{lastIt, curIt});
			}
			
			// Initiate first iteration
			CurrentIteration.Number = 1;

			while (true)
			{
				// Calculate element forces
				FemInput.Elements.CalculateForces();

				// Update residual
				ResidualUpdate();

				// Check convergence or if analysis must stop
				if (IterativeStop())
					return;
				
				// Add iteration
				_iterations.Add(new IterationResult(FemInput.NumberOfDoFs) { Number = CurrentIteration.Number});
				
				// Update stiffness and displacements
				UpdateStiffness();
				DisplacementUpdate();

				// Increase iteration count
				CurrentIteration.Number++;
			}
		}

		/// <summary>
		///     Check if iterative procedure must stop by achieving convergence or achieving the maximum number of iterations.
		/// </summary>
		/// <remarks>
		///     If the maximum number of iterations is reached, <see cref="Stop" /> is set to true.
		/// </remarks>
		/// <returns>
		///     True if convergence is reached or the maximum number of iterations is reached.
		/// </returns>
		private bool IterativeStop()
		{
			// Check if one stop condition is reached
			Stop = (int) CurrentIteration >= MaxIterations    || ResidualForces.ContainsNaNOrInfinity() ||
			       DisplacementVector!.ContainsNaNOrInfinity() || GlobalStiffness!.ContainsNaN();

			switch (Stop)
			{
				// Check if maximum number of iterations is reached
				case true:
					StopMessage = $"Convergence not reached at load step {(int) CurrentLoadStep}";
					return Stop;

				default:
					return
						VerifyConvergence(Convergence(CurrentIteration.ResidualForces, CurrentLoadStep.Forces));
			}
		}

		/// <summary>
		///     Update residual force <see cref="Vector{T}" />.
		/// </summary>
		private void ResidualUpdate() =>
			ResidualForces = FemInput.AssembleInternalForces() - CurrentLoadStep.Forces;

		/// <summary>
		///     Save load step results after achieving convergence.
		/// </summary>
		private void SaveLoadStepResults()
		{
			var curLoadStep = CurrentLoadStep;
			curLoadStep.IsCalculated  = true;
			curLoadStep.Displacements = DisplacementVector!;
			curLoadStep.Stiffness     = GlobalStiffness!;

			if (!_monitoredIndex.HasValue)
				return;

			// Get displacement
			var disp = Length.FromMillimeters(DisplacementVector![_monitoredIndex.Value]);

			// Set to load step
			curLoadStep.MonitoredDisplacement = new MonitoredDisplacement(disp, (double) (int) curLoadStep / NumLoadSteps);
		}

		/// <summary>
		///     Execute step by step analysis.
		/// </summary>
		/// <param name="simulate">Set true to execute analysis until convergence is not achieved (structural failure).</param>
		private void StepAnalysis(bool simulate)
		{
			// Initiate first load step
			var loadStep = 1;

			while (simulate || loadStep <= NumLoadSteps)
			{
				// Get the force vector
				var f = (double) loadStep / NumLoadSteps * ForceVector;

				// Create load step
				_loadSteps.Add(new LoadStepResult(loadStep, f));

				// Iterate
				Iterate();

				// Verify if convergence was not reached
				if (Stop)
				{
					CorrectResults();
					return;
				}

				// Set load step results
				SaveLoadStepResults();

				// Increment load step
				loadStep++;
			}
		}

		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="convergence">
		///     Calculated convergence.
		///     <para>See: <see cref="Convergence" />.</para>
		/// </param>
		private bool VerifyConvergence(double convergence) => convergence <= Tolerance && (int) CurrentIteration >= MinIterations;

		#endregion

	}
}