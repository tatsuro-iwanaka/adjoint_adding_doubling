#include<iostream>
#include<aad/App>

int main(int argc, char *argv[])
{
	if(argc != 2)
	{
		std::cout << "Usage: ./radiative_transfer <configuration filename>" << std::endl;

		return 1;
	}

	aad::RadiativeTransfer rtm(argv[1]);
	rtm.run();
	rtm.exportResultNetCDF();

	return 0;
}