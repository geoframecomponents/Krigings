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

import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;

import org.junit.Test;

import experimentalVariogram.ExperimentalVariogram;


public class TestExperimentalVariogram{
	@SuppressWarnings("nls")
	
	@Test
	public void testVariogram() throws Exception {


		//

		String stationIdField = "Id";

		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/experimentalVGM/jura.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file ="resources/Input/experimentalVGM/variogram_test.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";

		reader.initProcess();

		ExperimentalVariogram Meuse = new ExperimentalVariogram();

		Meuse.inStations = stationsFC;
		Meuse.fStationsid = stationIdField;


		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/experimentalVGM/experimental_distances.csv";
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		

		OmsTimeSeriesIteratorWriter writerS = new OmsTimeSeriesIteratorWriter();
		writerS.file = "resources/Output/experimentalVGM/experimental_variogram.csv";
		writerS.tStart = reader.tStart;
		writerS.tTimestep = reader.tTimestep;

		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			Meuse.inData = id2ValueMap;


			Meuse.process();

			HashMap<Integer, double[]> resultD = Meuse.outDistances;
			HashMap<Integer, double[]> resultS = Meuse.outExperimentalVariogram;


			writer.inData = resultD;
			writer.writeNextLine();

			writerS.inData = resultS;
			writerS.writeNextLine();
		}
		//
		reader.close();
		writer.close();
		writerS.close();


	}

}




