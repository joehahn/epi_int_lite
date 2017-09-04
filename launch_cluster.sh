#!/bin/bash 

#launch_cluster.sh
#by Joe Hahn, jmh.datasciences@gmail.com, 3 August 2017.
#
#launch a small EMR cluster 
#
#To execute:    ./launch_cluster.sh


#choose aws-cli profile (no dashes plz)
profile_str="spiralwaves"

#set aws region...note athena is available in us-west-2 (oregon) but not in us-west-1
aws_region="us-west-2"

#create s3 bucket
bucket_name="spiralwaves"
aws s3 mb "s3://$bucket_name" --profile "$profile_str"

#copy bootstrap script to s3
aws s3 cp bootstrap.sh "s3://$bucket_name/scripts/bootstrap.sh" --profile "$profile_str"

#copy piggyback script to s3
aws s3 cp piggyback.sh "s3://$bucket_name/scripts/piggyback.sh" --profile "$profile_str"

#copy provision_datasci script to s3
aws s3 cp provision_datasci.sh "s3://$bucket_name/scripts/provision_datasci.sh" --profile "$profile_str"

#copy aws access keys to s3, they will be needed by athena
aws s3 cp private/accessKeys.csv "s3://$bucket_name/accessKeys.csv" --profile "$profile_str"

#upload this repo to s3...I rather that the provision_datasci script clone this repo,
#but i couldn't resolve ssh issues...
cd ..
tar --exclude='spiral-waves/private' --exclude='spiral-waves/.git' \
    -zcvf /tmp/spiral-waves.tar.gz spiral-waves &> /dev/null
aws s3 cp /tmp/spiral-waves.tar.gz "s3://$bucket_name/spiral-waves.tar.gz" --profile "$profile_str"
cd spiral-waves

#select hadoop applications to be installed
#noting that spark also provides beeline that will be used to talk to athena
applications='Name=Hadoop Name=Spark'

#set ec2_attributes
ec2_attributes='{"KeyName":"datasci","InstanceProfile":"EMR_EC2_DefaultRole","SubnetId":"subnet-087c7641"}'

#set instance_groups
#m3.xlarge (4cpu & 15Gb) costs $6.38/day=$190/month
#m4.2xlarge (8cpu & 32Gb) costs $10.34/day=$310/month
#using EMR bumps cost up by 23%
#athena and s3 charges are negligible
instance_groups='[{"InstanceCount":1,"InstanceGroupType":"MASTER","InstanceType":"m4.2xlarge","Name":"Master - 1"},{"InstanceCount":2,"InstanceGroupType":"CORE","InstanceType":"m4.2xlarge","Name":"Core - 2"}]'

#launch the cluster...change to --no-auto-terminate to persist the cluster
echo 'launching EMR cluster...'
aws emr create-cluster \
    --profile $profile_str \
    --auto-scaling-role EMR_AutoScaling_DefaultRole \
    --applications $applications \
    --ec2-attributes "$ec2_attributes" \
    --service-role EMR_DefaultRole \
    --enable-debugging \
    --release-label emr-5.7.0 \
    --log-uri "s3n://$bucket_name/elasticmapreduce/" \
    --name "$profile_str" \
    --tags "Name=$profile_str" \
    --instance-groups "$instance_groups" \
    --region "$aws_region" \
    --bootstrap-action Path="s3://$bucket_name/scripts/bootstrap.sh" \
    --steps Type=CUSTOM_JAR,Name=CustomJAR,ActionOnFailure=CONTINUE,Jar=s3://$aws_region.elasticmapreduce/libs/script-runner/script-runner.jar,Args=["s3://$bucket_name/scripts/piggyback.sh"] \
    --no-auto-terminate \
    --no-termination-protected

##launch the datascience instance
#echo 'launching datasci instance...'
#./launch_datasci.sh
