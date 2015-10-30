#pragma once

int width(int edgeType)
{
	switch(edgeType)
	{
		case(1):			// motorway	
		case(2):
				return 10;
				break;
		case(3):			// primary
		case(4):
				return 6;
				break;
		case(5):			// secondary
		case(6):
				return 3;
				break;
		case(7):			// tertiary 
		case(8):
				return 2;
				break;
		case(9):			// trunk
		case(10):	return 8;
				break;

		default:
				return 1;
				break;
	}
}

int color(int edgeType)
{
	// z.B.			1=blau (Autobahnen)
	//			2=rot (Bundesstraßen)
	//			3=orange (Landstraßen)
	//			4=gelb (Kreisstraßen)
	//			5=schwarz (alle anderen)
	switch(edgeType)
	{
		case(1):			// motorway	
		case(2):
				return 1;
				break;
		case(3):			// primary
		case(4):
				return 2;
				break;
		case(5):			// secondary
		case(6):
				return 3;
				break;
		case(7):			// tertiary 
		case(8):
				return 4;
				break;

		case(9):			// trunk
		case(10):	return 2;
				break;
		default:
				return 5;
				break;
		
	}
}
