#include <ipcp_obj.h>

#define RINA_PREFIX "ipcm.mad.ribd_v1"
#include <librina/logs.h>
#include <librina/exceptions.h>

#include "../agent.h"
#include "../../../ipcm.h"
#include "../ribf.h"

namespace rinad {
namespace mad {
namespace rib_v1 {

//Static class names
const std::string IPCPObj::class_name = "IPCProcess";

//Class
IPCPObj::IPCPObj(int ipcp_id) : RIBObj(class_name), processID_(ipcp_id){

}

void IPCPObj::read(const rina::cdap_rib::con_handle_t &con,
				const std::string& fqn,
				const std::string& class_,
				const rina::cdap_rib::filt_info_t &filt,
				const int invoke_id,
				rina::cdap_rib::SerializedObject &obj_reply,
				rina::cdap_rib::res_info_t& res){

	(void) con;
	(void) fqn;
	(void) class_;
	(void) filt;
	(void) invoke_id;
	res.code_ = rina::cdap_rib::CDAP_SUCCESS;

	mad_manager::structures::ipcp_t info;
	info.process_id = processID_;
	info.name = IPCManager->get_ipcp_name(processID_);
	//TODO: Add missing stuff...

	encoder.encode(info, obj_reply);
}

bool IPCPObj::delete_(const rina::cdap_rib::con_handle_t &con,
				const std::string& fqn,
				const std::string& class_,
				const rina::cdap_rib::filt_info_t &filt,
				const int invoke_id,
				rina::cdap_rib::res_info_t& res){
	(void) con;
	(void) fqn;
	(void) class_;
	(void) filt;
	(void) invoke_id;

	//Fill in the response
	res.code_ = rina::cdap_rib::CDAP_SUCCESS;

	//Call the IPCManager and return
	if (IPCManager->destroy_ipcp(ManagementAgent::inst, processID_) != IPCM_SUCCESS) {
		LOG_ERR("Unable to destroy IPCP with id %d", processID_);
		res.code_ = rina::cdap_rib::CDAP_ERROR;
		return false;
	}

	return true;
}

void IPCPObj::create_cb(const rina::rib::rib_handle_t rib,
			const rina::cdap_rib::con_handle_t &con,
			const std::string& fqn,
			const std::string& class_,
			const rina::cdap_rib::filt_info_t &filt,
			const int invoke_id,
			const rina::cdap_rib::SerializedObject &obj_req,
			rina::cdap_rib::SerializedObject &obj_reply,
			rina::cdap_rib::res_info_t& res){

	IPCPObj* ipcp;

	(void)con;
	(void)filt;
	(void)invoke_id;
	(void)obj_reply;

	//Set return value
	res.code_ = rina::cdap_rib::CDAP_SUCCESS;

	if (class_ != IPCPObj::class_name) {
		LOG_ERR("Create operation failed: received an invalid class name '%s' during create operation in '%s'",
								class_.c_str(),
								fqn.c_str());
		res.code_ = rina::cdap_rib::CDAP_INVALID_OBJ_CLASS;
		return;
	}


	mad_manager::structures::ipcp_config_t object;
	rinad::mad_manager::encoders::IPCPConfigEncoder().decode(
			obj_req, object);

	//Call the IPCManager
	int ipcp_id = createIPCP(object);

	if (ipcp_id == -1){
		LOG_ERR("Create operation failed in '%s'. Unknown error.",
								fqn.c_str());

		res.code_ = rina::cdap_rib::CDAP_ERROR;
		return;
	}

	if(!assignToDIF(object, ipcp_id)){
		LOG_ERR("Create operation failed: unable to assign IPCP with id '%d' to DIF in path '%s'. Unknown error.",
								ipcp_id,
								fqn.c_str());
		res.code_ = rina::cdap_rib::CDAP_ERROR;
		return;
	}

	if (!object.dif_to_register.empty()) {
		if (!registerAtDIF(object, ipcp_id)){
			LOG_ERR("Create operation failed: unable to register IPCP with id '%d' in path '%s'. Unknown error.",
								ipcp_id,
								fqn.c_str());
			res.code_ = rina::cdap_rib::CDAP_ERROR;
			return;
		}
	}

	try {
		ipcp = new IPCPObj(ipcp_id);
	} catch (...) {
		LOG_ERR("Create operation failed for ipcp id '%d' in path '%s': out of memory.",
								ipcp_id,
								fqn.c_str());
		res.code_ = rina::cdap_rib::CDAP_ERROR;
		return;
	}

	//Finally add it into the RIB
	try{
		RIBFactory::getProxy()->addObjRIB(rib, fqn, &ipcp);
	}catch(...) {
		LOG_ERR("Create operation failed: for ipcp id '%d' in path '%s'. Out of memory.",
								ipcp_id,
								fqn.c_str());
		res.code_ = rina::cdap_rib::CDAP_ERROR;
		return;
	}
}

int IPCPObj::createIPCP(
		rinad::mad_manager::structures::ipcp_config_t &object) {

	CreateIPCPPromise ipcp_promise;

	rina::ApplicationProcessNamingInformation ipcp_name(
			object.process_name, object.process_instance);

	if (IPCManager->create_ipcp(ManagementAgent::inst, &ipcp_promise,
					ipcp_name, object.process_type)
			== IPCM_FAILURE
			|| ipcp_promise.wait() != IPCM_SUCCESS) {
		LOG_ERR("Error while creating IPC process");
		return -1;
	}
	LOG_INFO("IPC process created successfully [id = %d, name= %s]",
			ipcp_promise.ipcp_id, object.process_name.c_str());
	return ipcp_promise.ipcp_id;

}

bool IPCPObj::assignToDIF(
			rinad::mad_manager::structures::ipcp_config_t &object,
			int ipcp_id) {

	// ASSIGN TO DIF
	Promise assign_promise;

	rina::ApplicationProcessNamingInformation dif_name(
							object.dif_to_assign,
							std::string());

	if (!IPCManager->ipcp_exists(ipcp_id)) {
		LOG_ERR("No such IPC process id %d", ipcp_id);
		return false;
	}

	if (IPCManager->assign_to_dif(ManagementAgent::inst, &assign_promise,
					ipcp_id, dif_name) == IPCM_FAILURE
			|| assign_promise.wait() != IPCM_SUCCESS) {
		LOG_ERR("DIF assignment failed");
		return false;
	}

	LOG_INFO("DIF assignment completed successfully");

	return true;
}

bool IPCPObj::registerAtDIF(mad_manager::structures::ipcp_config_t &object,
								int ipcp_id) {

	Promise promise;

	rina::ApplicationProcessNamingInformation dif_name(
			object.dif_to_register, std::string());

	if (!IPCManager->ipcp_exists(ipcp_id)) {
		LOG_ERR("No such IPC process id");
		return false;
	}

	if (IPCManager->register_at_dif(ManagementAgent::inst, &promise,
					ipcp_id, dif_name) == IPCM_FAILURE
			|| promise.wait() != IPCM_SUCCESS) {
		LOG_ERR("Registration failed");
		return false;
	}

	LOG_INFO("IPC process registration completed successfully");

	return true;
}



};//namespace rib_v1
};//namespace mad
};//namespace rinad
