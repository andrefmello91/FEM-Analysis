using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Monitored displacement struct.
	/// </summary>
	public struct MonitoredDisplacement
	{
		/// <summary>
		///		Get the displacement value, at current <see cref="LoadFactor"/>.
		/// </summary>
		public Length Displacement { get; }
		
		/// <summary>
		///		Get the load factor associated to <see cref="Displacement"/>.
		/// </summary>
		public double LoadFactor { get; }

		/// <summary>
		///		Monitored displacement constructor.
		/// </summary>
		/// <param name="displacement">The displacement value, at current <paramref name="loadFactor"/>.</param>
		/// <param name="loadFactor">The load factor associated to <paramref name="displacement"/>.</param>
		public MonitoredDisplacement(Length displacement, double loadFactor)
		{
			Displacement = displacement;
			LoadFactor   = loadFactor;
		}

		public override string ToString() => 
			$"Displacement = {Displacement}\n" +
			$"Load Factor = {LoadFactor:0.00}";
	}
}