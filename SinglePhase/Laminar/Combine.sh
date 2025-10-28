for base in $(printf '%s\n' Dump_file/Output_field/cell_dump_t*_pid*.csv | sed -E 's/_pid[0-9]+\.csv$//' | sort -u); do
    awk 'FNR==1 && NR!=1 {next} 1' "${base}_pid"*.csv > "${base}_all.csv" && rm "${base}_pid"*.csv
done
