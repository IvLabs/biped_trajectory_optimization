function offset_origins = travel_offset(impacts, static_origins)

 
offset_origins = [static_origins, zeros(size(static_origins,1),2)];

for index = 1:size(impacts,1) 
    impact_index = impacts(index,2);

    offset_origins(impact_index+1:end,1:2:end) = ... 
        offset_origins(impact_index+1:end,1:2:end) ... 
        + static_origins(impact_index,7);

    offset_origins(impact_index+1:end,2:2:end) = ... 
        offset_origins(impact_index+1:end,2:2:end) ... 
        + static_origins(impact_index,8); 
  end

end

