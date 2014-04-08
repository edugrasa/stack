package rina.ipcprocess.impl.PDUForwardingTable.ribobjects;


import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import rina.cdap.api.CDAPSessionDescriptor;
import rina.cdap.api.message.CDAPMessage;
import rina.ipcprocess.api.IPCProcess;
import rina.ribdaemon.api.BaseRIBObject;
import rina.ribdaemon.api.ObjectInstanceGenerator;
import rina.ribdaemon.api.RIBDaemonException;
import rina.PDUForwardingTable.api.FlowStateObject;
import rina.PDUForwardingTable.api.FlowStateObjectGroup;
import rina.PDUForwardingTable.api.PDUFTable;

public class FlowStateRIBObjectGroup extends BaseRIBObject{
	private static final Log log = LogFactory.getLog(FlowStateRIBObjectGroup.class);
	
	private PDUFTable pdufTable = null;
	
	/*		Constructors		*/
	public FlowStateRIBObjectGroup(PDUFTable table, IPCProcess ipcProcess)
	{
		super(ipcProcess, FlowStateObjectGroup.FLOW_STATE_GROUP_RIB_OBJECT_CLASS, 
				ObjectInstanceGenerator.getObjectInstance(), FlowStateObjectGroup.FLOW_STATE_GROUP_RIB_OBJECT_NAME);
		pdufTable = table;
	}
	@Override
	public Object getObjectValue(){
		Object[] result = new Object[this.getChildren().size()];
		for(int i=0; i<result.length; i++){
			result[i] = this.getChildren().get(i).getObjectValue();
		}
		return result;
	}
	@Override
	public void write(CDAPMessage cdapMessage, CDAPSessionDescriptor cdapSessionDescriptor)
			throws RIBDaemonException 
	{
		pdufTable.writeMessageRecieved(cdapMessage, cdapSessionDescriptor.getPortId());
	}
	
	@Override
	public void create(String objectClass, long objectInstance, String objectName, Object object) throws RIBDaemonException{
		if (!(object instanceof FlowStateObject)){
			throw new RIBDaemonException(RIBDaemonException.OBJECTCLASS_DOES_NOT_MATCH_OBJECTNAME, 
					"Object class ("+object.getClass().getName()+") does not match object name "+objectName);
		}
		
		FlowStateRIBObject ribObject = new FlowStateRIBObject((FlowStateObject) object, getIPCProcess(), objectName);
		this.addChild(ribObject);
		getRIBDaemon().addRIBObject(ribObject);
	}
	
	@Override
	public void delete(Object object) throws RIBDaemonException 
	{
		FlowStateObject fso = (FlowStateObject)object;
		this.removeChild(fso.getID());
		getRIBDaemon().removeRIBObject(fso.getID());
	}
	

}