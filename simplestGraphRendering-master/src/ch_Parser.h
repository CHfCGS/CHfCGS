#pragma once

#include "nodesAndEdges.h"

class CH_Parser{
		
	void createCHNode(std::string input_string, CHNode& r_node) {
		std::stringstream ss(input_string);
		std::string buffer;
		
		ss >> buffer; //fmiID
		ss >> buffer; //osmID
		
		ss >> buffer;
		r_node.lat = atof(buffer.c_str());

		ss >> buffer;
		r_node.lon = atof(buffer.c_str());
		
		ss >> buffer; //elev
		
		ss >> buffer;
		r_node.lvl = atof(buffer.c_str());			
	}
	
	void createCHEdge(std::string input_string, CHEdge& r_edge) {
		std::stringstream ss(input_string);
		std::string buffer;
		
		ss >> buffer;
		r_edge.src = atoi(buffer.c_str());
		
		ss >> buffer;
		r_edge.tgt = atoi(buffer.c_str());
		
		ss >> buffer;
		//r_edge.dist = atoi(buffer.c_str());
		
		ss >> buffer;
		r_edge.type = atoi(buffer.c_str());
		
		ss >> buffer;
		r_edge.speed = atoi(buffer.c_str());
						
		ss >> buffer;
		r_edge.child_edge1 = atof(buffer.c_str());

		ss >> buffer;
		r_edge.child_edge2 = atof(buffer.c_str());
		/*
		is >> src;
        is >> tgt;
        is >> dist;
        is >> type;
        is >> speed;
        is >> child_edge1;                
        is >> child_edge2;      
		*/
		/*
		inputFile>>curNodeExt.fmiID;
				inputFile>>curNodeExt.osmID;
				inputFile>>curNode.lat;
				inputFile>>curNode.lon;
				inputFile>>curNodeExt.elev;
				
				inputFile>>srcID;
				// curEdge.source=srcID;	// change here depending on desired space efficiency
				inputFile>>curEdge.target;
				inputFile>>curEdge.weight;	
				if (curEdge.weight==0)
				{
					countZeroWeight++;
				}
				inputFile>>curEdge.type;
				inputFile.getline(junkC,256);
				* */
	}
	
	public:
		bool parseTxtGraphFile(std::string graphfile, std::vector<CHNode>& n, std::vector<CHEdge>& e) {
			std::string buffer;
			std::ifstream file;

			file.open(graphfile.c_str(), std::ios::in);

			if( file.is_open())
			{
				file.seekg(0, std::ios::beg);
				
				//junk
				for (int i = 0; i < 10; i++){
					getline(file,buffer,'\n');     
				}
				
				getline(file,buffer,'\n');
				uint node_count = (uint)atoi(buffer.c_str());
				getline(file,buffer,'\n');
				uint edge_count = (uint)atoi(buffer.c_str());

				n.reserve(node_count);
				for(uint i=0; i<node_count; i++)
				{
					getline(file,buffer,'\n');
					n.push_back(CHNode());
					createCHNode(buffer, n.back() );
				}

				e.reserve(edge_count);
				for(uint j=0; j<edge_count; j++)
				{
					getline(file,buffer,'\n');
					e.push_back(CHEdge());
					createCHEdge(buffer, e.back() );
				}
				file.close();

				return true;
			}

			return false;
		}
};
