function [Channel] = CreateChannel( NoSources, NoDestinations, AntennaPerSource, AntennaPerDestination)

% This function is to generate the channel
% NoSources : Number of source objects
% NoDestinations: Number of destination objects
% AntennaPerSource: Number of antennas per source (scalar or vector)
% AntennaPerDestination: Number of antennas per destination (scalar or vector)


if (length(AntennaPerSource)>1)
	if (length(AntennaPerSource)<NoSources)
		disp('Error: Length of vector providing the number of anntennas per source is not enough');
		exit(0)
	end
else
	AntennaPerSource = repmat(AntennaPerSource,1,NoSources);
end

if (length(AntennaPerDestination)>1)
	if (length(AntennaPerDestination)<NoDestinations)
		disp('Error: Length of vector providing the number of anntennas per destination is not enough');
		exit(0)
	end
else
	AntennaPerDestination = repmat(AntennaPerDestination,1,NoDestinations);
end

for i = 1:1:NoSources
	for j=1:1:NoDestinations

		Channel(:,:,i,j) = sqrt(1/2)*(randn(AntennaPerSource(i),AntennaPerDestination(j)) + 1i*randn(AntennaPerSource(i),AntennaPerDestination(j)) );

	end
end
	
	
end