/* This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
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
package variogram;

import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;

import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.modules.ModelsEngine;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;

import org.jgrasstools.hortonmachine.i18n.HortonMessageHandler;
import org.opengis.feature.simple.SimpleFeature;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

@Description("Experimental semivariogram algorithm.")
@Documentation("Experimental semivariogram")
@Author(name = "Giuseppe Formetta, Francesco Adami, Silvia Franceschi & Marialaura Bancheri")
@Keywords("Experimental semivariogram, Kriging, Hydrology")
@Label(JGTConstants.STATISTICS)
@Name("variogram")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class Variogram extends JGTModel {

	@Description("The vector of the measurement point, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;

	@Description("The field of the vector of stations, defining the elevation.")
	@In
	public String fStationsZ = null;

	@Description("The file with the measured data, to be interpolated.")
	@In
	public HashMap<Integer, double[]> inData = null;

	@Description("The path to the printed file.")
	@In
	public String pPath = null;

	@Description("Spatial separation distance up to which point pairs are included in semivariance estimates; "
			+ "as a default, the length of the diagonal of the box spanning the data is divided by three.")
	@In
	public double pCutoff;

	@Description("The Experimental Variogram.")
	@Out
	public double[][] outResult = null;

	@Description("The Experimental Distances.")
	@Out
	public double[] outDist = null;

	@Description("The Experimental Variogram.")
	@Out
	public double[] outVar = null;

	@Description("The Experimental Pairs Number.")
	@Out
	public double[] outNumPairs = null;

	@Description("All are equal.")
	@Out
	public boolean outAllEquals = true;


	@Description("All are different.")
	@Out
	public int outAllDiff = 0;

	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new LogProgressMonitor();

	private HortonMessageHandler msg = HortonMessageHandler.getInstance();

	public boolean areAllEquals;

	public int differents;

	@Execute
	public void process() throws Exception {

		List<Double> xStationList = new ArrayList<Double>();
		List<Double> yStationList = new ArrayList<Double>();
		List<Double> zStationList = new ArrayList<Double>();
		List<Double> hStationList = new ArrayList<Double>();

		/*
		 * counter for the number of station with measured value !=0.
		 */
		int n1 = 0;

		/*
		 * Store the station coordinates and measured data in the array.
		 */
		FeatureIterator<SimpleFeature> stationsIter = inStations.features();
		try {
			while (stationsIter.hasNext()) {
				SimpleFeature feature = stationsIter.next();
				int id = ((Number) feature.getAttribute(fStationsid))
						.intValue();
				double z = 0;
				if (fStationsZ != null) {
					try {
						z = ((Number) feature.getAttribute(fStationsZ))
								.doubleValue();
					} catch (NullPointerException e) {
						pm.errorMessage(msg.message("kriging.noStationZ"));
						throw new Exception(msg.message("kriging.noStationZ"));

					}
				}
				Coordinate coordinate = ((Geometry) feature
						.getDefaultGeometry()).getCentroid().getCoordinate();


				// h is the vector with measured data
				double[] h = inData.get(id);
				if (h == null || isNovalue(h[0])) {

					/*
					 * skip data for non existing stations.
					 * Also skip novalues.
					 */
					continue;
				}

				if (Math.abs(h[0]) >= 0.0) { // TOLL
					xStationList.add(coordinate.x);
					yStationList.add(coordinate.y);
					zStationList.add(z);
					hStationList.add(h[0]);
					n1 = n1 + 1;
				}

			}
		} finally {
			stationsIter.close();
		}

		int nStaz = xStationList.size();


		/*
		 * The coordinates of the station points plus in the last position a place
		 * for the coordinate of the point to interpolate.
		 */
		double[] xStation = new double[nStaz];
		double[] yStation = new double[nStaz];
		double[] zStation = new double[nStaz];
		double[] hStation = new double[nStaz];

		areAllEquals = true;

		differents = 0;

		if (nStaz != 0) {
			xStation[0] = xStationList.get(0);
			yStation[0] = yStationList.get(0);
			zStation[0] = zStationList.get(0);
			hStation[0] = hStationList.get(0);

			double previousValue = hStation[0];

			for (int i = 1; i < nStaz; i++) {

				double xTmp = xStationList.get(i);
				double yTmp = yStationList.get(i);
				double zTmp = zStationList.get(i);
				double hTmp = hStationList.get(i);

				boolean doubleStation = ModelsEngine.verifyDoubleStation(
						xStation, yStation, zStation, hStation, xTmp, yTmp,
						zTmp, hTmp, i, false, pm);

				if (!doubleStation) {
					xStation[i] = xTmp;
					yStation[i] = yTmp;
					zStation[i] = zTmp;
					hStation[i] = hTmp;
					if (areAllEquals && hStation[i] != previousValue) {

						areAllEquals = false;
					}
					if (hStation[i] != previousValue) {
						differents += 1;
					}
					previousValue = hStation[i];
				}
			}
		}


		if (differents > 2) {

			outResult = processAlgorithm(xStation, yStation, hStation, pCutoff);

			if (pPath != null && pPath.length() > 0) {
				FileWriter Rstatfile = new FileWriter(pPath);
				PrintWriter errestat = new PrintWriter(Rstatfile);
				for (int i = 0; i < (outResult.length + 1); i++) {

					for (int j = 0; j < (outResult[0].length); j++) {
						if (i == 0) {
							errestat.print("Np" + " " + "Dist" + " " + "Gamma"
									+ " " + "Moran" + " " + "Geary");
							break;
						}

						errestat.print(outResult[i - 1][j] + " ");

					}
					errestat.println();
					System.out.println();
				}

				Rstatfile.close();
			}
		}
		else{
			System.out.println("Only 1 data >0 or All the data are equal. Variogram is not running");
		}
		outAllEquals = areAllEquals;
		outAllDiff = differents;

	}

	public double[][] processAlgorithm(double[] xcord, double ycoord[],
			double[] values, double Cutoffinput) {

		double x1, x2, y1, y2;
		double dDifX, dDifY;
		double value;
		double mean = 0;
		double maxDistance = 0;

		int Cutoff_divide = 15;
		double Cutoff;
		int iCount = xcord.length;
		double distanceMatrix[][] = new double[iCount][iCount];
		double x_max = xcord[0], y_max = ycoord[0], diagonal;
		double x_min = xcord[0], y_min = ycoord[0];
		for (int i = 1; i < iCount; i++) {

			x_min = Math.min(x_min, xcord[i]);
			y_min = Math.min(y_min, ycoord[i]);
			x_max = Math.max(x_max, xcord[i]);
			y_max = Math.max(y_max, ycoord[i]);

		}

		diagonal = Math.sqrt((x_max - x_min) * (x_max - x_min)
				+ (y_max - y_min) * (y_max - y_min));

		if (Cutoffinput == 0) {
			Cutoff = diagonal / 3;
		} else
			Cutoff = Cutoffinput;


		// Compute the distance matrix
		for (int i = 0; i < iCount; i++) {
			x1 = xcord[i];
			y1 = ycoord[i];
			value = values[i];

			mean += value;

			for (int j = 0; j < iCount; j++) {

				x2 = xcord[j];
				y2 = ycoord[j];

				dDifX = x2 - x1;
				dDifY = y2 - y1;

				// Pitagora theorem
				distanceMatrix[i][j] = Math.sqrt(dDifX * dDifX + dDifY * dDifY); 

				maxDistance = Math.max(maxDistance, distanceMatrix[i][j]);

			}
		}

		// compute the mean of the input values
		mean /= (double) iCount; 

		double[][] result = calculate(Cutoff_divide, Cutoff, distanceMatrix, values, mean, maxDistance);

		return result;

	} 

	public double[][] calculate(int num, double cutoff, double[][] distanceMatrix, double[] values, double mean,
			double maxDistance) {

		int iClasses;
		int iClass;

		double binAmplitude = cutoff / num;


		// number of distance for each bin
		iClasses = (int) (maxDistance / binAmplitude + 2); 

		// definition of the vectors containing the variance, covariance, semivariance, 
		// 	number of the points in the specified bin..
		double[] m_dMoran = new double[iClasses];

		double[] m_dGeary = new double[iClasses];

		double[] dDen = new double[iClasses];

		double[] m_dSemivar = new double[iClasses]; 

		int[] iPointsInClass = new int[iClasses]; 

		boolean bIsInClass[] = new boolean[iClasses];

		double[] m_ddist = new double[iClasses];


		for (int i = 0; i < distanceMatrix.length; i++) {
			Arrays.fill(bIsInClass, false); 

			// first cycle input values 
			double value1 = values[i]; 

			for (int j = i + 1; j < distanceMatrix.length; j++) {
				if (distanceMatrix[i][j] > 0 && distanceMatrix[i][j] < cutoff) {

					// return the class of considered distance 
					iClass = (int) Math.floor((distanceMatrix[i][j]) / binAmplitude); 

					// counts the number of distances of each class
					iPointsInClass[iClass]++; 

					// second cycle input values
					double value2 = values[j]; 

					// compute the numerator of the semivariance 
					double dSemivar = Math.pow((value1 - value2), 2.); 

					// sum all the semivariances for the considered class
					m_dSemivar[iClass] += dSemivar; 

					// compute the numerator of the Moran, 
					// summing all the covariances for the considered class
					m_dMoran[iClass] += (value1 - mean) * (value2 - mean);

					// compute the numerator of the Geary
					m_dGeary[iClass] = m_dSemivar[iClass]; 

					// flag the considered class as true since it has been processed 
					bIsInClass[iClass] = true; 

					m_ddist[iClass] += distanceMatrix[i][j]; 

				}
			}

			for (int j = 0; j < iClasses; j++) {
				if (bIsInClass[j]) { 

					// compute  the sum of all the variances for the class
					dDen[j] += Math.pow(value1 - mean, 2.); 

				}
			}

		}

		double[][] result = new double[iClasses][5];
		int countNONzero = 0;	
		int[] indicesNONzero = new int[iClasses];

		for (int u = 0; u < dDen.length; u++) {
			if (dDen[u] > 0) {
				countNONzero += 1;
				indicesNONzero[countNONzero - 1] = u;
			}
		}

		for (int i = 0; i < iClasses; i++) {
			if (dDen[i] != 0) {

				// Compute the Moran for each class
				m_dMoran[i] /= dDen[i];

				// Compute the Geary for each class
				m_dGeary[i] *= ((iPointsInClass[i] - 1) / (2. * iPointsInClass[i] * dDen[i])); // n-1/2n(var)

				// Compute the semivariance
				m_dSemivar[i] /= (2. * iPointsInClass[i]); 

				// Compute the mean distance for each class
				m_ddist[i] /= iPointsInClass[i];

				result[i][0] = iPointsInClass[i];
				result[i][1] = m_ddist[i];
				result[i][2] = m_dSemivar[i];
				result[i][3] = m_dMoran[i];
				result[i][4] = m_dGeary[i];

			}
		}
		outDist = new double[countNONzero];
		outVar = new double[countNONzero];
		double results[][] = new double[countNONzero][5];

		outNumPairs = new double[countNONzero];
		if (areAllEquals == false && differents > 2) {

			for (int ii = 0; ii < countNONzero; ii++) {

				results[ii][0] = result[indicesNONzero[ii]][0];
				outNumPairs[ii] = result[indicesNONzero[ii]][0];

				results[ii][1] = result[indicesNONzero[ii]][1];
				outDist[ii] = result[indicesNONzero[ii]][1];

				results[ii][2] = result[indicesNONzero[ii]][2];
				outVar[ii] = result[indicesNONzero[ii]][2];

				results[ii][3] = result[indicesNONzero[ii]][3];
				results[ii][4] = result[indicesNONzero[ii]][4];
			}
		} else {
			results = new double[iClasses][5];
			outNumPairs = new double[iClasses];
			outDist = new double[iClasses];
			outVar = new double[iClasses];

			for (int ii = 0; ii < iClasses; ii++) {
				results[ii][0] = JGTConstants.doubleNovalue;
				outNumPairs[ii] = JGTConstants.doubleNovalue;

				results[ii][1] = JGTConstants.doubleNovalue;
				outDist[ii] = JGTConstants.doubleNovalue;

				results[ii][2] = JGTConstants.doubleNovalue;
				outVar[ii] = JGTConstants.doubleNovalue;

				results[ii][3] = JGTConstants.doubleNovalue;
				results[ii][4] = JGTConstants.doubleNovalue;
			}

		}
		return results;

	}

}