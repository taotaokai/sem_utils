
job_id=0
for work in green misfit srcfrechet dxs dmt search
do 

  job_file=$event_dir/$work.job


  sbatch -j after_ok:$jid $job_file > $log_sbatch
  grep -i error $event_dir/green.sub > /dev/null
  if [ $? -eq 0 ]
  then
    echo "[ERROR] failed to submit $event_dir/green.job"
    cat $event_dir/green.sub
  fi

  job_id=$(grep "Submitted batch job" $event_dir/green.sub | awk '{print $NF}')

done
