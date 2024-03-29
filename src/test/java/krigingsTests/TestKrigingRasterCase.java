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
package krigingsTests;

import java.util.HashMap;


import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.hortonmachine.gears.io.rasterreader.OmsRasterReader;
import org.hortonmachine.gears.io.rasterwriter.OmsRasterWriter;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;

import org.junit.Test;
import org.junit.Assert;



import krigingsRasterCase.Krigings;
import pathGenerator.PathGenerator;



public class TestKrigingRasterCase {


	/**
	 * Run the kriging models.
	 *
	 * <p>
	 * This is the case which all the station have the same value.
	 * </p>
	 * @throws Exception
	 * @throws Exception
	 */

	@Test 
	public void testKriging() throws Exception {
		//
		String stationIdField = "ID_P";

		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/krigings/RasterCase/final_station.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		// OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		// interpolatedPointsReader.file = puntiFile.getAbsolutePath();
		// interpolatedPointsReader.readFeatureCollection();

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = "resources/Input/krigings/RasterCase/P1994_2013.csv";
		reader.idfield = "ID_P";
		reader.tStart = "1994-01-01 00:00";
		reader.tTimestep = 60;
		 reader.tEnd = "1994-01-01 01:00";
		reader.fileNovalue = "-999";
		reader.initProcess();

		Krigings kriging = new Krigings();
		
		OmsRasterReader demReader = new OmsRasterReader();
		demReader.file = "resources/Input/krigings/RasterCase/dem.asc";

		demReader.process();
		GridCoverage2D dem = demReader.outRaster;

		
		
		kriging.inGridCoverage2D = dem;

		kriging.inStations = stationsFC;
		kriging.fStationsid = stationIdField;

		kriging.range = 123537.0;
		kriging.nugget = 0.0;
		kriging.sill= 1.678383;
		kriging.pSemivariogramType="linear";
		
		PathGenerator path=new PathGenerator();


		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			/*
			 * Extract the result.
			 */

			path.pathToOutData="resources/Output/krigings/RasterCase/krigings.asc";
			path.tCurrent=reader.tCurrent;
			path.process();
			


			GridCoverage2D krigingRaster = kriging.outGrid;
			
			OmsRasterWriter writerRaster = new OmsRasterWriter();
			writerRaster.inRaster = krigingRaster;
			writerRaster.file = path.pathOutDataComplete;
			writerRaster.process();

			
		}

		reader.close();
		// writer.close();
	}
}
