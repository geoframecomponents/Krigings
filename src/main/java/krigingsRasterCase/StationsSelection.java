/*
 * GNU GPL v3 License
 *
 * Copyright 2016 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package krigingsRasterCase;

import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.jgrasstools.gears.libs.modules.ModelsEngine;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.gears.utils.sorting.QuickSortAlgorithm;
import org.opengis.feature.simple.SimpleFeature;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

// TODO: Auto-generated Javadoc
/**
 * The Class StationsSelection.
 */
public class StationsSelection {

	/** The .shp of the measurement point, containing the position of the stations. */
	public SimpleFeatureCollection inStations = null;

	/** The HM with the measured data to be interpolated. */
	public HashMap<Integer, double[]> inData = null;


	/** Include zeros in computations (default is true). */
	public boolean doIncludezero = true;

	/** In the case of kriging with neighbor, maxdist is the maximum distance 
	    within the algorithm has to consider the stations */
	public double maxdist;

	/** In the case of kriging with neighbor, inNumCloserStations is the number 
    of stations the algorithm has to consider */
	public int inNumCloserStations;

	/** The field of the vector of stations, defining the id.*/
	public String fStationsid = null;


	/** The pm. */
	public IJGTProgressMonitor pm = new LogProgressMonitor();

	/** The x station initial set. */
	public double[] xStationInitialSet;

	/** The y station initial set. */
	public double[] yStationInitialSet;

	/** The z station initial set. */
	public double[] zStationInitialSet;

	/** The h station initial set. */
	public double[] hStationInitialSet;

	/** The id station initial set. */
	public int[] idStationInitialSet;
	
	/** The flag is true if all the stations are equals */
	public boolean areAllEquals = true;
	
	/** The number of different stations. */
	public int n1 = 0;
	
	
	/** The model selectio for the choice of the stations. */
	Model modelSelection;




	/**
	 * Execute.
	 *
	 * @throws Exception the exception
	 */
	public void execute() throws Exception {

		// create the arraylist containing the station with the measurements 
		List<Double> xStationList = new ArrayList<Double>();
		List<Double> yStationList = new ArrayList<Double>();
		List<Double> zStationList = new ArrayList<Double>();
		List<Double> hStationList = new ArrayList<Double>();
		List<Integer> idStationList = new ArrayList<Integer>();


		/*
		 * Store the station coordinates and measured data in the array.
		 * Skip data for non existing stations and also skip novalues.
		 */

		FeatureIterator<SimpleFeature> stationsIter = inStations.features();
		try {
			while (stationsIter.hasNext()) {
				SimpleFeature feature = stationsIter.next();
				int id = ((Number) feature.getAttribute(fStationsid)).intValue();

				double z = 0;

				Coordinate coordinate = ((Geometry) feature.getDefaultGeometry()).getCentroid().getCoordinate();
				double[] h = inData.get(id);
				if (h == null || isNovalue(h[0])) {

					continue;
				}


				if (doIncludezero) {
					if (Math.abs(h[0]) >= 0.0) { // TOLL
						xStationList.add(coordinate.x);
						yStationList.add(coordinate.y);
						zStationList.add(z);
						hStationList.add(h[0]);
						idStationList.add(id);
						n1 = n1 + 1;
					}

				} else {
					if (Math.abs(h[0]) > 0.0) { // TOLL
						xStationList.add(coordinate.x);
						yStationList.add(coordinate.y);
						zStationList.add(z);
						hStationList.add(h[0]);
						idStationList.add(id);

						n1 = n1 + 1;
					}
				}
			}
		} finally {
			stationsIter.close();
		}

		int nStaz = xStationList.size();


		/*
		 * Check if the coordinates or the values are the same for all the measurements stations.
		 * xStationInitialSet has the dimensions of the coordinates of the measurements points 
		 * plus 1 (the station where it is going to interpolate)
		 */


		xStationInitialSet = new double[nStaz+ 1];
		yStationInitialSet = new double[nStaz+ 1];
		zStationInitialSet = new double[nStaz+ 1];
		hStationInitialSet = new double[nStaz+ 1];
		idStationInitialSet = new int[nStaz + 1];


		if (nStaz != 0) {
			xStationInitialSet[0] = xStationList.get(0);
			yStationInitialSet[0] = yStationList.get(0);
			zStationInitialSet[0] = zStationList.get(0);
			hStationInitialSet[0] = hStationList.get(0);
			idStationInitialSet[0] = idStationList.get(0);
			double previousValue = hStationInitialSet[0];

			/* for each station added to the vector, it checks the coordinate/values and if they are
			 * not null or equal, it adds to the list of the availble stations. If the coordinates
			 * or the values are all different, the flag areAllEquals  becomes false.   
			 */

			for (int i = 0; i < nStaz; i++) {

				double xTmp = xStationList.get(i);
				double yTmp = yStationList.get(i);
				double zTmp = zStationList.get(i);
				double hTmp = hStationList.get(i);
				int idTmp = idStationList.get(i);

				boolean doubleStation = ModelsEngine.verifyDoubleStation( xStationInitialSet, yStationInitialSet, zStationInitialSet, 
						hStationInitialSet,xTmp,yTmp, zTmp, hTmp, i, false, pm);
				if (!doubleStation) {
					xStationInitialSet[i] = xTmp;
					yStationInitialSet[i] = yTmp;
					zStationInitialSet[i] = zTmp;
					hStationInitialSet[i] = hTmp;
					idStationInitialSet[i] = idTmp;
					if (areAllEquals && hStationInitialSet[i] != previousValue) {
						areAllEquals = false;
					}
					if (hStationInitialSet[i] != previousValue) {
					}
					previousValue = hStationInitialSet[i];
				}
			}
		}

		/* in case of kriging with neighbor computes the distances between the
		 * point where is going to interpolate and the other stations and it
		 * sorts them
		 */ 

		if (inNumCloserStations > 0 || maxdist>0) {

			inNumCloserStations= (inNumCloserStations> nStaz)? nStaz:inNumCloserStations;	

			double x2, y2;
			double dDifX, dDifY;
			double distanceVector[] = new double[xStationInitialSet.length];
			double pos[] = new double[xStationInitialSet.length];


			for (int jj = 0; jj < xStationInitialSet.length; jj++) {

				x2 = xStationInitialSet[jj];
				y2 = yStationInitialSet[jj];

				dDifX = xStationInitialSet[n1] - x2;
				dDifY = yStationInitialSet[n1] - y2;
				distanceVector[jj] = Math.sqrt(dDifX * dDifX + dDifY * dDifY); 
				pos[jj] = jj;					
			}

			// sorts the distances
			QuickSortAlgorithm t = new QuickSortAlgorithm(pm);
			t.sort(distanceVector, pos);
			
			
			/*
			 * The dimension of the new vector of the station is then defined
			 * by the actual number of the station within the distance or defined 
			 * by the users
			 */

		
			modelSelection=SimpleModelFactory.createModel(distanceVector, inNumCloserStations, maxdist);
			int dim=modelSelection.numberOfStations();


			double[] xStationWithNeighbour = new double[dim];
			double[] yStationWithNeighbour = new double[dim];
			int[] idStationWithNeighbour = new int[dim];
			double[] zStationWithNeighbour = new double[dim];
			double[] hWithNeighbour = new double[dim];


			// it is necessary to actualize the counter of the stations
			n1=0;

			for (int i = 1; i < dim; i++) {					
				if (doIncludezero) {
					if (Math.abs(hStationInitialSet[(int) pos[i]]) >= 0.0) { // TOLL
						xStationWithNeighbour[n1] = xStationInitialSet[(int) pos[i]];
						yStationWithNeighbour[n1] = yStationInitialSet[(int) pos[i]];
						zStationWithNeighbour[n1] = zStationInitialSet[(int) pos[i]];
						idStationWithNeighbour[n1] = idStationInitialSet[(int) pos[i]];

						hWithNeighbour[n1] = hStationInitialSet[(int) pos[i]];
						n1 += 1;
					}
				} else {
					if (Math.abs(hStationInitialSet[(int) pos[i]]) > 0.0) {
						xStationWithNeighbour[n1] = xStationInitialSet[(int) pos[i]];
						yStationWithNeighbour[n1] = yStationInitialSet[(int) pos[i]];
						zStationWithNeighbour[n1] = zStationInitialSet[(int) pos[i]];
						idStationWithNeighbour[n1] = idStationInitialSet[(int) pos[i]];

						hWithNeighbour[n1] = hStationInitialSet[(int) pos[i]];
						n1 += 1;

					}
				}

			}

			xStationInitialSet = xStationWithNeighbour;
			yStationInitialSet = yStationWithNeighbour;
			zStationInitialSet = zStationWithNeighbour;
			idStationInitialSet = idStationWithNeighbour;
			hStationInitialSet = hWithNeighbour;

		}

	}

}