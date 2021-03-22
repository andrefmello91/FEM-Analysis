using System;
using System.Collections.Generic;
using System.Linq;
using Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet;
using UnitsNet.Units;
using static Extensions.UnitExtensions;

#nullable disable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Linear analysis class.
	/// </summary>
	public class Analysis
	{
		#region Fields

		/// <summary>
		///     Get the <see cref="FEMAnalysis.InputData" />.
		/// </summary>
		public InputData InputData { get; }

		#endregion

		#region Properties

		/// <inheritdoc cref="InputData.ConstraintIndex" />
		public int[] ConstraintIndex => InputData.ConstraintIndex;

		/// <summary>
		///     Get/set the displacement <see cref="Vector" />.
		/// </summary>
		/// <remarks>
		///     Components in <see cref="LengthUnit.Millimeter" />.
		/// </remarks>
		public Vector<double> DisplacementVector { get; protected set; }

		/// <inheritdoc cref="InputData.ForceVector" />
		public Vector<double> ForceVector { get; protected set; }

		/// <summary>
		///     Get/set global stiffness <see cref="Matrix" />.
		/// </summary>
		public Matrix<double> GlobalStiffness { get; protected set; }

		/// <summary>
		///     Get maximum force in <see cref="Stringers" /> and <see cref="Panels" />.
		/// </summary>
		public Force MaxElementForce => Max(MaxStringerForce, MaxPanelForce);

		/// <summary>
		///     Get maximum <see cref="Panels" /> force.
		/// </summary>
		public Force MaxPanelForce => Panels.Length > 0 ? Panels.Select(pnl => pnl.MaxForce).Max() : Force.Zero;

		/// <summary>
		///     Get maximum <see cref="Stringers" /> force.
		/// </summary>
		public Force MaxStringerForce => Stringers.Length > 0 ? Stringers.Select(str => str.MaxForce).Max() : Force.Zero;

		/// <summary>
		///     Get nodes of SPM model.
		/// </summary>
		public Node[] Nodes => InputData.Grips;

		/// <inheritdoc cref="InputData.NumberOfDoFs" />
		public int NumberOfDoFs => InputData.NumberOfDoFs;

		/// <summary>
		///     Get panels of SPM model.
		/// </summary>
		public Panel[] Panels => InputData.Panels;

		/// <summary>
		///     Get stringers of SPM model.
		/// </summary>
		public Stringer[] Stringers => InputData.Elements;

		#endregion

		#region Constructors

		/// <summary>
		///     Analysis base object.
		/// </summary>
		/// <param name="inputData">The <see cref="FEMAnalysis.InputData" /> for SPM analysis.</param>
		public Analysis(InputData inputData) => InputData = inputData;

		#endregion

		#region  Methods



		/// <summary>
		///     Calculate displacement <see cref="Vector" />.
		/// </summary>
		/// <param name="globalStiffness">Current global stiffness <see cref="Matrix" />.</param>
		/// <param name="forceVector">Current force <see cref="Vector" />.</param>
		public static Vector<double> CalculateDisplacements(Matrix<double> globalStiffness, Vector<double> forceVector) => globalStiffness.Solve(forceVector);

		/// <summary>
		///     Add element internal forces to global force <see cref="Vector" />.
		/// </summary>
		/// <param name="element">The <see cref="IFiniteElement" /> to add to <paramref name="globalForceVector" />.</param>
		/// <param name="globalForceVector">The global force <see cref="Vector" />.</param>
		public static void AddForce(IFiniteElement element, Vector<double> globalForceVector)
		{
			if (element is null)
				return;

			var dofIndex = element.DoFIndex;

			for (var i = 0; i < dofIndex.Length; i++)
			{
				// DoF index
				var j = dofIndex[i];

				// Add values
				globalForceVector[j] += element.Forces[i];
			}
		}

		/// <param name="elements">The <see cref="IFiniteElement" />'s to add to <paramref name="globalForceVector" />.</param>
		/// <inheritdoc cref="AddForce(IFiniteElement, Vector{double})" />
		public static void AddForce(IEnumerable<IFiniteElement> elements, Vector<double> globalForceVector)
		{
			if (elements is null || !elements.Any())
				return;

			foreach (var element in elements)
				AddForce(element, globalForceVector);
		}

		/// <summary>
		///     Do the analysis.
		/// </summary>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis.ForceVector" /> (default: 1).</param>
		public void Do(double loadFactor = 1)
		{
			// Set force vector
			ForceVector = InputData.ForceVector * loadFactor;

			// Assemble and simplify global stiffness and force vector
			UpdateStiffness();

			// Solve
			DisplacementVector = CalculateDisplacements(GlobalStiffness, ForceVector);

			// Calculate element displacements and forces
			ElementAnalysis(DisplacementVector);

			// Set nodal displacements
			NodalDisplacements(DisplacementVector);
		}

		/// <summary>
		///     Calculate force vector including reactions.
		/// </summary>
		public Vector<double> ForcesAndReactions()
		{
			// Calculate forces and reactions
			var f = GlobalStiffness * DisplacementVector;
			f.CoerceZero(1E-9);

			return f;
		}

		/// <summary>
		///     Update <see cref="GlobalStiffness" />.
		/// </summary>
		/// <param name="simplify">Simplify stiffness and force vector? (default: true)</param>
		protected void UpdateStiffness(bool simplify = true)
		{
			// Initialize the global stiffness matrix
			GlobalStiffness = AssembleStiffness();

			// Simplify stiffness matrix
			if (simplify)
				SimplifyStiffnessMatrix();
		}

		/// <summary>
		///     Assemble global stiffness <see cref="Matrix" />.
		/// </summary>
		protected Matrix<double> AssembleStiffness()
		{
			var stiffness = Matrix<double>.Build.Dense(NumberOfDoFs, NumberOfDoFs);

			AddStiffness(Stringers, stiffness);

			AddStiffness(Panels, stiffness);

			return stiffness;
		}

		/// <summary>
		///     Simplify <see cref="GlobalStiffness" /> and <see cref="ForceVector" />.
		/// </summary>
		protected void SimplifyStiffnessMatrix()
		{
			foreach (var i in ConstraintIndex)
			{
				// Clear the row and column [i] in the stiffness matrix (all elements will be zero)
				GlobalStiffness.ClearRow(i);
				GlobalStiffness.ClearColumn(i);

				// Set the diagonal element to 1
				GlobalStiffness[i, i] = 1;

				// Clear the row in the force vector
				if (ForceVector != null)
					ForceVector[i] = 0;

				// So ui = 0
			}

			// Simplification for internal nodes (There is only a displacement at the Stringer direction, the perpendicular one will be zero)
			foreach (var node in Nodes)
			{
				if (node.Type != NodeType.Internal)
					continue;

				// Get DoF indexes
				var index = node.DoFIndex;

				// Verify rows
				foreach (var i in index)
				{
					// Verify what line of the matrix is composed of zeroes
					if (GlobalStiffness.Row(i).Exists(num => !num.ApproxZero()))
						continue;

					// The row is composed of only zeroes, so the displacement must be zero
					// Set the diagonal element to 1
					GlobalStiffness[i, i] = 1;

					// Clear the row in the force vector
					if (ForceVector != null)
						ForceVector[i] = 0;
				}
			}

			// Approximate small numbers to zero
			GlobalStiffness.CoerceZero(1E-9);
		}

		/// <summary>
		///     Set displacements to <see cref="Stringers" /> and <see cref="Panels" /> and calculate element forces.
		/// </summary>
		/// <param name="globalDisplacements">
		///     Current displacement <see cref="Vector" />, with components in
		///     <see cref="LengthUnit.Millimeter" />.
		/// </param>
		protected void ElementAnalysis(Vector<double> globalDisplacements)
		{
			foreach (var stringer in Stringers)
				stringer.Analysis(globalDisplacements);

			foreach (var panel in Panels)
				panel.Analysis(globalDisplacements);
		}

		/// <summary>
		///     Set displacements to <see cref="Nodes" />.
		/// </summary>
		/// <param name="globalDisplacements">Current displacement <see cref="Vector" />.</param>
		protected void NodalDisplacements(Vector<double> globalDisplacements)
		{
			foreach (var node in Nodes)
				node.SetDisplacements(globalDisplacements);
		}

		/// <summary>
		///     Assemble the internal force <see cref="Vector" />.
		/// </summary>
		protected Vector<double> InternalForces()
		{
			var fi = Vector<double>.Build.Dense(NumberOfDoFs);

			AddForce(Stringers, fi);

			AddForce(Panels, fi);

			// Simplify for constraints
			foreach (var i in ConstraintIndex)
				fi[i] = 0;

			return fi;
		}


		/// <summary>
		///     Get indexes of continued stringers.
		/// </summary>
		protected List<(int str1, int str2)> ContinuedStringers()
		{
			// Initialize a Tuple to store the continued stringers
			var contStrs = new List<(int str1, int str2)>();

			// Calculate the parameter of continuity
			var par = 0.5 * Math.Sqrt(2);

			// Verify in the list what stringers are continuous
			foreach (var str1 in Stringers)
			foreach (var str2 in Stringers)
			{
				// Verify if it's other Stringer
				if (str1.Number == str2.Number)
					continue;

				// Create a tuple with the minimum Stringer number first
				var contStr = (Math.Min(str1.Number, str2.Number), Math.Max(str1.Number, str2.Number));

				// Verify if it's already on the list
				if (contStrs.Contains(contStr))
					continue;

				// Verify the cases
				// Case 1: stringers initiate or end at the same node
				if (str1.Grips[0] == str2.Grips[0] || str1.Grips[2] == str2.Grips[2])
				{
					// Get the direction cosines
					var (l1, m1) = str1.Geometry.Angle.DirectionCosines();
					var (l2, m2) = str2.Geometry.Angle.DirectionCosines();

					// Calculate the condition of continuity
					var cont = l1 * l2 + m1 * m2;

					// Verify the condition
					if (cont < -par) // continued Stringer
						contStrs.Add(contStr);
				}

				// Case 2: a Stringer initiate and the other end at the same node
				else if (str1.Grips[0] == str2.Grips[2] || str1.Grips[2] == str2.Grips[0])
				{
					// Get the direction cosines
					var (l1, m1) = str1.Geometry.Angle.DirectionCosines();
					var (l2, m2) = str2.Geometry.Angle.DirectionCosines();

					// Calculate the condition of continuity
					var cont = l1 * l2 + m1 * m2;

					// Verify the condition
					if (cont > par) // continued Stringer
						contStrs.Add(contStr);
				}
			}

			// Order the list
			contStrs = contStrs.OrderBy(str => str.Item2).ThenBy(str => str.Item1).ToList();

			// Return the list
			return contStrs;
		}

		/// <summary>
		///     Get the list of panel's DoFs that have continuity.
		/// </summary>
		private List<int> PanelContinuousDoFs()
		{
			var contDofs = new List<int>();

			foreach (var pnl1 in Panels)
			{
				// Get Dofs normal to edges
				var DoFs1 = pnl1.DoFIndex;
				var nDofs1 = new[]
				{
					DoFs1[1], DoFs1[2], DoFs1[5], DoFs1[6]
				};

				foreach (var pnl2 in Panels)
					if (pnl1.Number != pnl2.Number)
					{
						// Get Dofs normal to edges
						var DoFs2 = pnl2.DoFIndex;
						var nDofs2 = new[]
						{
							DoFs2[1], DoFs2[2], DoFs2[5], DoFs2[6]
						}.ToList();

						// Verify if there is a common DoF
						foreach (var dof in nDofs1)
							if (nDofs2.Contains(dof) && !contDofs.Contains(dof))
								contDofs.Add(dof);
					}
			}

			return
				contDofs;
		}

		public override string ToString() =>
			$"{InputData}\n" +
			"Global Stifness:\n" +
			$"{GlobalStiffness}\n" +
			"Displacement Vector:\n" +
			$"{DisplacementVector}";

		#endregion
	}
}